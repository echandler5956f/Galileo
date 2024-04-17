#include "galileo/legged-model/LeggedInterface.h"

namespace
{
    Eigen::VectorXd ReadVector(std::string str_vector)
    {
        // Remove the leading and trailing double quotes
        int idx_end = str_vector.find_last_of(')');
        int idx_start = str_vector.find_last_of('(');

        str_vector = str_vector.substr(idx_start + 1, idx_end - idx_start - 1);
        // Split the string by commas
        std::vector<std::string> str_values;
        std::stringstream ss(str_vector);
        std::string str_value;
        while (std::getline(ss, str_value, ','))
        {
            str_values.push_back(str_value);
        }

        // Convert the string values to double and store in Eigen::VectorXd
        Eigen::VectorXd vector(str_values.size());
        for (size_t i = 0; i < str_values.size(); ++i)
        {
            vector(i) = std::stod(str_values[i]);
        }

        return vector;
    }
}

namespace galileo
{
    namespace legged
    {
        LeggedInterface::LeggedInterface(std::string sol_data_dir, std::string plot_dir)
        {
            surfaces_ = std::make_shared<galileo::legged::environment::EnvironmentSurfaces>();
            solution_interface_ = std::make_shared<galileo::opt::solution::Solution>();
            meshcat_interface = std::make_shared<galileo::tools::MeshcatInterface>(sol_data_dir);
            plotting_interface = std::make_shared<galileo::tools::GNUPlotInterface>(plot_dir);
        }

        void LeggedInterface::LoadModel(std::string model_file_location, std::vector<std::string> end_effector_names)
        {
            std::string *end_effector_names_array = new std::string[end_effector_names.size()];

            for (size_t i = 0; i < end_effector_names.size(); ++i)
            {
                end_effector_names_array[i] = end_effector_names[i];
            }
            // std::vector<casadi::Dict> legged_options = {legged_opts_, legged_opts_, legged_opts_, legged_opts_};
            std::vector<casadi::Dict> legged_options = {casadi::Dict(), casadi::Dict(), casadi::Dict(), casadi::Dict()};

            robot_ = std::make_shared<LeggedBody>(model_file_location, end_effector_names.size(), end_effector_names_array, legged_options);
            states_ = robot_->si;
            model_file_location_ = model_file_location;
        }

        void LeggedInterface::LoadParameters(std::string parameter_file_location)
        {
            // Load the parameters from the given parameter file.
            std::map<std::string, std::tuple<std::string, std::string>> imported_vars;
            std::ifstream file(parameter_file_location);
            // Assert that file exists
            assert(file.is_open());
            std::string row;

            while (std::getline(file, row))
            {
                std::stringstream ss(row);
                std::string key;
                std::string value;
                std::string type;
                std::getline(ss, key, '|');
                std::getline(ss, value, '|');
                std::getline(ss, type, '\n');
                imported_vars[key] = std::make_tuple(value, type);
            }

            for (auto &var : imported_vars)
            {
                // If the first word up until "." is "constraints" then it is a constraint parameter for the legged model
                if (var.first.substr(0, var.first.find(".")) == "constraints")
                {
                    // Remove the "constraints." prefix
                    std::string key = var.first.substr(var.first.find(".") + 1);
                    // Remove the "constraints." prefix
                    std::string value = std::get<0>(var.second);
                    // Remove the "constraints." prefix
                    std::string type = std::get<1>(var.second);
                    if (type == "double")
                    {
                        constraint_params_[key] = std::stod(value);
                        std::cout << "Constraint parameter: " << key << " = " << constraint_params_[key] << std::endl;
                    }
                }

                // If the first word up until "." is "nlp" then it is a nlp parameter
                if (var.first.substr(0, var.first.find(".")) == "nlp")
                {
                    // Remove the "nlp." prefix
                    std::string key = var.first.substr(var.first.find(".") + 1);
                    // Remove the "nlp." prefix
                    std::string value = std::get<0>(var.second);
                    // Remove the "nlp." prefix
                    std::string type = std::get<1>(var.second);
                    if (type == "int")
                        opts_[key] = std::stoi(value);
                    else if (type == "double")
                        opts_[key] = std::stod(value);
                    else if (type == "string")
                        opts_[key] = value;
                    else if (type == "bool")
                        opts_[key] = (value == "true");
                    else if (type == "vector")
                        continue;
                }
            }

            // Extract the string value from imported_vars
            Eigen::VectorXd Q_diag = ReadVector(std::get<0>(imported_vars["cost.Q_diag"]));
            Eigen::VectorXd R_diag = ReadVector(std::get<0>(imported_vars["cost.R_diag"]));

            casadi::DM joint_lb = casadi::DM::zeros(states_->nvju, 1);
            casadi::DM joint_ub = casadi::DM::zeros(states_->nvju, 1);
            Eigen::VectorXd joint_lb_vec = robot_->model.lowerPositionLimit.block(states_->nqb, 0, states_->nvju, 1);
            Eigen::VectorXd joint_ub_vec = robot_->model.upperPositionLimit.block(states_->nqb, 0, states_->nvju, 1);
            tools::eigenToCasadi(joint_lb_vec, joint_lb);
            tools::eigenToCasadi(joint_ub_vec, joint_ub);

            assert(joint_lb.rows() == states_->nvju);
            assert(joint_ub.rows() == states_->nvju);

            joint_limits_.lower = joint_lb;
            joint_limits_.upper = joint_ub;

            assert(Q_diag.size() == states_->ndx);
            assert(R_diag.size() == states_->nF + states_->nF);

            cost_params_.Q_diag = Q_diag;
            cost_params_.R_diag = R_diag;

            if (imported_vars.find("cost.terminal_weight") != imported_vars.end())
            {
                cost_params_.terminal_weight = std::stod(std::get<0>(imported_vars["cost.terminal_weight"]));
                std::cout << "Terminal weight: " << cost_params_.terminal_weight << std::endl;
            }

            if (imported_vars.find("solver") != imported_vars.end())
                solver_type_ = std::get<0>(imported_vars["solver"]);

            parameters_set_ = true;
        }

        void LeggedInterface::CreateProblemData(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            assert(robot_ != nullptr); // Model must be loaded
            casadi::Function Phi;

            CreateCost(initial_state, target_state, Phi);

            std::shared_ptr<opt::GeneralProblemData> gp_data = std::make_shared<opt::GeneralProblemData>(robot_->fint, robot_->fdiff, Phi);

            problem_data_ = std::make_shared<LeggedRobotProblemData>(gp_data,
                                                                     surfaces_,
                                                                     robot_->contact_sequence,
                                                                     states_, std::make_shared<legged::ADModel>(robot_->cmodel),
                                                                     std::make_shared<legged::ADData>(robot_->cdata),
                                                                     robot_->getEndEffectors(),
                                                                     robot_->cx, robot_->cu, robot_->cdt, initial_state, target_state, joint_limits_, constraint_params_);
        }

        // TODO: Generate a reference trajectory somewhere and share it betwen the objective and initial guess
        void LeggedInterface::CreateCost(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state, casadi::Function &Phi)
        {
            assert(robot_ != nullptr);

            casadi::DM X0 = initial_state;

            assert(cost_params_.Q_diag.size() == states_->ndx);
            Eigen::MatrixXd Q_mat = cost_params_.Q_diag.asDiagonal();

            assert(cost_params_.R_diag.size() == states_->nF + states_->nF);
            Eigen::MatrixXd R_taskspace = cost_params_.R_diag.asDiagonal();

            casadi::DM q0_dm = states_->get_q(X0);
            Eigen::VectorXd q0 = Eigen::Map<Eigen::VectorXd>(q0_dm.get_elements().data(), q0_dm.size1() * q0_dm.size2());

            Eigen::MatrixXd R_mat = robot_->initializeInputCostWeight(R_taskspace, q0);

            casadi::SX Q = casadi::SX::zeros(states_->ndx, states_->ndx);
            casadi::SX R = casadi::SX::zeros(states_->nu, states_->nu);

            pinocchio::casadi::copy(Q_mat, Q);
            pinocchio::casadi::copy(R_mat, R);

            casadi::SX Xf_sym = casadi::SX::sym("Xf_sym", states_->nx, 1);

            /*Both f_state_error and fdiff work for the quaternion error, but fdiff uses the quaternion logarithm, so it is more accurate albeit slightly more expensive*/
            casadi::SX X_error = robot_->f_state_error(casadi::SXVector{robot_->cx, Xf_sym}).at(0);
            // casadi::SX X_error = robot_->fdiff(casadi::SXVector{rrobot_->cx, Xf_sym, 1.}).at(0);

            pinocchio::computeTotalMass(robot_->cmodel, robot_->cdata);

            for (std::size_t i = 0; i < robot_->contact_sequence->getPhases().size(); ++i)
            {
                casadi::SX U_ref = robot_->weightCompensatingInputsForPhase(i);
                casadi::SX u_error = robot_->cu - U_ref;
                casadi::Function L = casadi::Function("L_" + std::to_string(i),
                                                      {robot_->cx, robot_->cu, Xf_sym},
                                                      {0.5 * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error)) +
                                                       0.5 * casadi::SX::dot(u_error, casadi::SX::mtimes(R, u_error))});
                robot_->contact_sequence->FillPhaseCost(i, L);
            }

            // TODO: Add the terminal cost weight to a parameter file
            Phi = casadi::Function("Phi",
                                   {robot_->cx, Xf_sym},
                                   {cost_params_.terminal_weight * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error))});
        }

        std::vector<LeggedInterface::LeggedConstraintBuilderType>
        LeggedInterface::getLeggedConstraintBuilders() const
        {
            LeggedConstraintBuilderType friction_cone_constraint_builder =
                std::make_shared<constraints::FrictionConeConstraintBuilder<LeggedRobotProblemData>>();

            LeggedConstraintBuilderType velocity_constraint_builder =
                std::make_shared<constraints::VelocityConstraintBuilder<LeggedRobotProblemData>>();

            LeggedConstraintBuilderType contact_constraint_builder =
                std::make_shared<constraints::ContactConstraintBuilder<LeggedRobotProblemData>>();

            return {velocity_constraint_builder, friction_cone_constraint_builder};
        }

        void LeggedInterface::setContactSequence(std::vector<int> knot_num, std::vector<double> knot_time, std::vector<uint> mask_vec, std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces)
        {
            assert(robot_ != nullptr); // Model must be loaded
            std::shared_ptr<contact::ContactSequence> contact_sequence = std::make_shared<contact::ContactSequence>(robot_->num_end_effectors_);

            for (size_t j = 0; j < mask_vec.size(); ++j)
            {
                contact::ContactMode mode;
                mode.combination_definition = robot_->getContactCombination(mask_vec[j]);
                mode.contact_surfaces = contact_surfaces[j];
                contact_sequence->addPhase(mode, knot_num[j], knot_time[j]);
            }

            setContactSequence(contact_sequence);
        }

        void LeggedInterface::Initialize(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            assert(CanInitialize());

            // Create the problem data from the loaded parameter values
            CreateProblemData(initial_state, target_state);

            // Create the trajectory optimizer.
            CreateTrajOpt();

            fully_initialized_ = true;
        }

        void LeggedInterface::CreateTrajOpt()
        {
            auto constraint_builders = getLeggedConstraintBuilders();
            decision_builder_ = std::make_shared<galileo::legged::constraints::LeggedDecisionDataBuilder<LeggedRobotProblemData>>();

            std::lock_guard<std::mutex> lock(trajectory_opt_mutex_);

            trajectory_opt_ = std::make_shared<LeggedTrajOpt>(problem_data_, robot_->contact_sequence, constraint_builders, decision_builder_, opts_, solver_type_);

            trajectory_opt_->InitFiniteElements(3);
        }

        void LeggedInterface::Update(double global_time, const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            if (!first_update_)
            {
                std::lock_guard<std::mutex> lock_sol(solution_mutex_);
                GetNextGuess(global_time);
            }

            // Solve the problem
            {
                std::lock_guard<std::mutex> lock_traj(trajectory_opt_mutex_);
                if (first_update_)
                {
                    trajectory_opt_->AdvanceFiniteElements(global_time, initial_state, target_state);
                }
                else
                {
                    casadi::DM w0;
                    tools::eigenToCasadi(curr_guess_, w0);
                    trajectory_opt_->AdvanceFiniteElements(global_time, initial_state, target_state, w0);
                }
            }

            {
                std::lock_guard<std::mutex> lock_sol(solution_mutex_);
                //@todo Akshay5312, reevaluate thread safety
                solution_interface_->UpdateSolution(trajectory_opt_->getSolutionSegments());
                solution_interface_->UpdateConstraints(trajectory_opt_->getConstraintDataSegments());
            }

            curr_time_ = global_time;
        }

        void LeggedInterface::GetNextGuess(double global_time)
        {
            double dt = global_time - curr_time_;
            casadi::DMVector stacked_decision_variable_times = trajectory_opt_->getStackedDecisionVariableTimes();
            solution_interface_->GetNextGuessFromPrevSolution(stacked_decision_variable_times, curr_guess_, dt);
        }

        // void LeggedInterface::UpdateProblemBoundaries(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        // {
        //     problem_data_->legged_decision_problem_data.X0 = (initial_state);

        //     // Create a new cost
        //     casadi::Function Phi;
        //     CreateCost(initial_state, target_state, Phi);

        //     problem_data_->gp_data->Phi = Phi;
        // }

        bool LeggedInterface::GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result)
        {
            std::lock_guard<std::mutex> lock_sol(solution_mutex_);
            return solution_interface_->GetSolution(query_times, state_result, input_result);
        }

        void LeggedInterface::PlotTrajectories(const Eigen::VectorXd &query_times, const Eigen::MatrixXd &state_result, const Eigen::MatrixXd &input_result)
        {
            opt::solution::solution_t solution(query_times, state_result.transpose(), input_result.transpose());

            // Collect the data specific to each end effector
            std::vector<std::tuple<size_t, size_t>> wrench_indices;
            std::vector<std::vector<std::string>> wrench_legend_names;
            std::vector<std::string> ee_plot_names;
            for (auto ee : robot_->getEndEffectors())
            {
                wrench_indices.push_back(states_->frame_id_to_index_range[ee.second->frame_id]);
                if (ee.second->is_6d)
                    wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}", "\\tau_{x}", "\\tau_{y}", "\\tau_{z}"});
                else
                {
                    wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}"});
                }
                ee_plot_names.push_back("Contact Wrench of " + ee.second->frame_name);
            }

            plotting_interface->PlotSolution(solution, {std::make_tuple(states_->q_index, states_->q_index + 3), std::make_tuple(states_->q_index + 3, states_->q_index + states_->nqb)},
                                             {wrench_indices},
                                             {"Positions", "Orientations"},
                                             {{"x", "y", "z"}, {"qx", "qy", "qz", "qw"}},
                                             {ee_plot_names},
                                             {wrench_legend_names});
        }

        void LeggedInterface::VisualizeSolutionAndConstraints(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result)
        {
            std::unique_lock<std::mutex> lock_sol(solution_mutex_);
            std::vector<std::vector<galileo::opt::constraint_evaluations_t>> constraints = solution_interface_->GetConstraints(query_times, state_result, input_result);
            lock_sol.unlock();

            Eigen::MatrixXd subMatrix = state_result.block(states_->q_index, 0, states_->nq, state_result.cols());

            meshcat_interface->WriteTimes(query_times, "sol_times.csv");
            meshcat_interface->WriteJointPositions(subMatrix, "sol_states.csv");
            meshcat_interface->WriteMetadata(model_file_location_, "metadata.csv");

            opt::solution::solution_t solution(query_times, state_result.transpose(), input_result.transpose());

            // Collect the data specific to each end effector
            std::vector<std::tuple<size_t, size_t>> wrench_indices;
            std::vector<std::vector<std::string>> wrench_legend_names;
            std::vector<std::string> ee_plot_names;
            for (auto ee : robot_->getEndEffectors())
            {
                wrench_indices.push_back(states_->frame_id_to_index_range[ee.second->frame_id]);
                if (ee.second->is_6d)
                    wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}", "\\tau_{x}", "\\tau_{y}", "\\tau_{z}"});
                else
                {
                    wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}"});
                }
                ee_plot_names.push_back("Contact Wrench of " + ee.second->frame_name);
            }

            plotting_interface->PlotSolution(solution, {std::make_tuple(states_->q_index, states_->q_index + 3), std::make_tuple(states_->q_index + 3, states_->q_index + states_->nqb)},
                                             {wrench_indices},
                                             {"Positions", "Orientations"},
                                             {{"x", "y", "z"}, {"qx", "qy", "qz", "qw"}},
                                             {ee_plot_names},
                                             {wrench_legend_names});

            plotting_interface->PlotConstraints(constraints);
        }
    }
}
