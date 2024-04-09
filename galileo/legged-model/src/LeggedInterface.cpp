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
            // TODO: Add these options to a parameter file
            casadi::Dict legged_opts;
            // // legged_opts["cse"] = true;
            // legged_opts["jit"] = true;
            // legged_opts["jit_options.flags"] = "-Ofast -march=native -ffast-math";
            // legged_opts["jit_options.compiler"] = "gcc";
            // // legged_opts["jit_options.temp_suffix"] = false;
            // legged_opts["compiler"] = "shell";
            // // legged_opts["jit_cleanup"] = false;

            robot_ = std::make_shared<LeggedBody>(model_file_location, end_effector_names.size(), end_effector_names_array, legged_opts);
            states_ = robot_->si;
            model_file_location_ = model_file_location;
        }

        void LeggedInterface::LoadParameters(std::string parameter_file_location)
        {
            // Load the parameters from the given parameter file.
            std::map<std::string, std::string> imported_vars;
            std::ifstream file(parameter_file_location);
            // assert that file exists
            assert(file.is_open());
            std::string row;

            while (std::getline(file, row))
            {
                std::stringstream ss(row);
                std::string key;
                std::string value;
                std::getline(ss, key, ',');
                std::getline(ss, value, '\n');
                imported_vars[key] = value;
            }

            /*Extract the constraint parameters. TODO: Find a better way to do this...*/

            if (imported_vars.find("mu") != imported_vars.end())
                constraint_params_["mu"] = std::stod(imported_vars["mu"]);
            if (imported_vars.find("normal_force_max") != imported_vars.end())
                constraint_params_["normal_force_max"] = std::stod(imported_vars["normal_force_max"]);
            if (imported_vars.find("ideal_offset_height") != imported_vars.end())
                constraint_params_["ideal_offset_height"] = std::stod(imported_vars["ideal_offset_height"]);
            if (imported_vars.find("footstep_height_scaling") != imported_vars.end())
                constraint_params_["footstep_height_scaling"] = std::stod(imported_vars["footstep_height_scaling"]);
            if (imported_vars.find("max_following_leeway_planar") != imported_vars.end())
                constraint_params_["max_following_leeway_planar"] = std::stod(imported_vars["max_following_leeway_planar"]);
            if (imported_vars.find("min_following_leeway_planar") != imported_vars.end())
                constraint_params_["min_following_leeway_planar"] = std::stod(imported_vars["min_following_leeway_planar"]);
            if (imported_vars.find("footstep_vel_start") != imported_vars.end())
                constraint_params_["footstep_vel_start"] = std::stod(imported_vars["footstep_vel_start"]);
            if (imported_vars.find("footstep_vel_end") != imported_vars.end())
                constraint_params_["footstep_vel_end"] = std::stod(imported_vars["footstep_vel_end"]);

            /*Extract the solver parameters. TODO: Find a better way to do this...*/

            solver_type_ = imported_vars["solver"];

            if (solver_type_ == "ipopt")
            {
                if (imported_vars.find("ipopt.fixed_variable_treatment") != imported_vars.end())
                    opts_["ipopt.fixed_variable_treatment"] = imported_vars["ipopt.fixed_variable_treatment"];
                if (imported_vars.find("ipopt.max_iter") != imported_vars.end())
                    opts_["ipopt.max_iter"] = std::stoi(imported_vars["ipopt.max_iter"]);
                if (imported_vars.find("ipopt.linear_solver") != imported_vars.end())
                    opts_["ipopt.linear_solver"] = imported_vars["ipopt.linear_solver"];
                if (imported_vars.find("ipopt.hessian_approximation") != imported_vars.end())
                    opts_["ipopt.hessian_approximation"] = imported_vars["ipopt.hessian_approximation"];
                if (imported_vars.find("ipopt.ma97_order") != imported_vars.end())
                    opts_["ipopt.ma97_order"] = imported_vars["ipopt.ma97_order"];
                if (imported_vars.find("pass_nonlinear_variables") != imported_vars.end())
                    opts_["pass_nonlinear_variables"] = (imported_vars["pass_nonlinear_variables"] == "true");

                if (imported_vars.find("ipopt.limited_memory_aug_solver ") != imported_vars.end())
                    opts_["ipopt.limited_memory_aug_solver"] = imported_vars["ipopt.limited_memory_aug_solver"];
                if (imported_vars.find("ipopt.limited_memory_max_history") != imported_vars.end())
                    opts_["ipopt.limited_memory_max_history"] = std::stoi(imported_vars["ipopt.limited_memory_max_history"]);
                if (imported_vars.find("ipopt.limited_memory_update_type") != imported_vars.end())
                    opts_["ipopt.limited_memory_update_type"] = imported_vars["ipopt.limited_memory_update_type"];
                if (imported_vars.find("ipopt.limited_memory_initialization") != imported_vars.end())
                    opts_["ipopt.limited_memory_initialization"] = imported_vars["ipopt.limited_memory_initialization"];
                if (imported_vars.find("ipopt.limited_memory_init_val") != imported_vars.end())
                    opts_["ipopt.limited_memory_init_val"] = std::stod(imported_vars["ipopt.limited_memory_init_val"]);

                if (imported_vars.find("ipopt.limited_memory_init_val_max") != imported_vars.end())
                    opts_["ipopt.limited_memory_init_val_max"] = std::stod(imported_vars["ipopt.limited_memory_init_val_max"]);
                if (imported_vars.find("ipopt.limited_memory_init_val_min") != imported_vars.end())
                    opts_["ipopt.limited_memory_init_val_min"] = std::stod(imported_vars["ipopt.limited_memory_init_val_min"]);
                if (imported_vars.find("ipopt.limited_memory_max_skipping") != imported_vars.end())
                    opts_["ipopt.limited_memory_max_skipping"] = std::stoi(imported_vars["ipopt.limited_memory_max_skipping"]);
                if (imported_vars.find("ipopt.limited_memory_special_for_resto") != imported_vars.end())
                    opts_["ipopt.limited_memory_special_for_resto"] = imported_vars["ipopt.limited_memory_special_for_resto"];
                if (imported_vars.find("ipopt.hessian_approximation_space") != imported_vars.end())
                    opts_["ipopt.hessian_approximation_space"] = imported_vars["ipopt.hessian_approximation_space"];
            }
            else if (solver_type_ == "snopt")
            {
                opts_["snopt.Major iterations limit"] = std::stoi(imported_vars["snopt.Major iterations limit"]);
                opts_["snopt.Minor iterations limit"] = std::stoi(imported_vars["snopt.Minor iterations limit"]);
                opts_["snopt.Iterations limit"] = std::stoi(imported_vars["snopt.Iterations limit"]);
                opts_["snopt.Cold Start/Warm Start"] = imported_vars["snopt.Cold Start/Warm Start"];
            }

            // Extract the string value from imported_vars
            Eigen::VectorXd Q_diag = ReadVector(imported_vars["Q_diag"]);
            Eigen::VectorXd R_diag = ReadVector(imported_vars["R_diag"]);

            assert(Q_diag.size() == states_->ndx);
            assert(R_diag.size() == states_->nu);

            cost_params_.Q_diag = Q_diag;
            cost_params_.R_diag = R_diag;
            if (imported_vars.find("terminal_weight") != imported_vars.end())
                cost_params_.terminal_weight = std::stod(imported_vars["terminal_weight"]);

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
                                                                     robot_->cx, robot_->cu, robot_->cdt, initial_state, target_state, constraint_params_);
        }

        // TODO: Generate a reference trajectory somewhere and share it betwen the objective and initial guess
        void LeggedInterface::CreateCost(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state, casadi::Function &Phi)
        {
            assert(robot_ != nullptr);

            casadi::DM X0 = initial_state;

            assert(cost_params_.Q_diag.size() == states_->ndx);
            Eigen::MatrixXd Q_mat = cost_params_.Q_diag.asDiagonal();

            assert(cost_params_.R_diag.size() == states_->nu);
            Eigen::MatrixXd R_taskspace = cost_params_.R_diag.asDiagonal();

            casadi::DM q0_dm = states_->get_q(X0);
            Eigen::VectorXd q0 = Eigen::Map<Eigen::VectorXd>(q0_dm.get_elements().data(), q0_dm.size1() * q0_dm.size2());

            Eigen::MatrixXd R_mat = robot_->initializeInputCostWeight(R_taskspace, q0);

            casadi::SX Q = casadi::SX::zeros(states_->ndx, states_->ndx);
            casadi::SX R = casadi::SX::zeros(states_->nu, states_->nu);

            pinocchio::casadi::copy(Q_mat, Q);
            pinocchio::casadi::copy(R_mat, R);

            /*Both f_state_error and fdiff work for the quaternion error, but fdiff uses the quaternion logarithm, so it is more accurate albeit slightly more expensive*/
            casadi::SX X_error = robot_->f_state_error(casadi::SXVector{robot_->cx, target_state}).at(0);
            // casadi::SX X_error = robot_->fdiff(casadi::SXVector{rrobot_->cx, target_state, 1.}).at(0);

            pinocchio::computeTotalMass(robot_->cmodel, robot_->cdata);

            for (std::size_t i = 0; i < robot_->contact_sequence->getPhases().size(); ++i)
            {
                casadi::SX U_ref = robot_->weightCompensatingInputsForPhase(i);
                casadi::SX u_error = robot_->cu - U_ref;
                casadi::Function L = casadi::Function("L_" + std::to_string(i),
                                                      {robot_->cx, robot_->cu},
                                                      {0.5 * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error)) +
                                                       0.5 * casadi::SX::dot(u_error, casadi::SX::mtimes(R, u_error))});
                robot_->contact_sequence->FillPhaseCost(i, L);
            }

            // TODO: Add the terminal cost weight to a parameter file
            Phi = casadi::Function("Phi",
                                   {robot_->cx},
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
        }

        void LeggedInterface::Update(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            UpdateProblemBoundaries(initial_state, target_state);

            // Solve the problem
            std::lock_guard<std::mutex> lock_traj(trajectory_opt_mutex_);
            trajectory_opt_->optimize();

            std::lock_guard<std::mutex> lock_sol(solution_mutex_);
            //@todo Akshay5312, reevaluate thread safety
            solution_interface_->UpdateSolution(trajectory_opt_->getSolutionSegments());

            solution_interface_->UpdateConstraints(trajectory_opt_->getConstraintDataSegments());
        }

        bool LeggedInterface::GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result)
        {
            std::lock_guard<std::mutex> lock_sol(solution_mutex_);
            return solution_interface_->GetSolution(query_times, state_result, input_result);
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

        void LeggedInterface::UpdateProblemBoundaries(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            problem_data_->legged_decision_problem_data.X0 = (initial_state);

            // Create a new cost
            casadi::Function Phi;
            CreateCost(initial_state, target_state, Phi);

            problem_data_->gp_data->Phi = Phi;

            std::lock_guard<std::mutex> lock(trajectory_opt_mutex_);
            trajectory_opt_->initFiniteElements(1, initial_state);
        }

    }
}
