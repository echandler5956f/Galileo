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
        LeggedInterface::LeggedInterface()
        {
            surfaces_ = std::make_shared<galileo::legged::environment::EnvironmentSurfaces>();
            solution_interface_ = std::make_shared<galileo::opt::solution::Solution>();
            meshcat_interface = std::make_shared<galileo::tools::MeshcatInterface>("../examples/visualization/solution_data/");
            plotting_interface = std::make_shared<galileo::tools::GNUPlotInterface>("../examples/visualization/plots/");
        }

        void LeggedInterface::LoadModel(std::string model_file_location, std::vector<std::string> end_effector_names)
        {
            std::string *end_effector_names_array = new std::string[end_effector_names.size()];

            for (size_t i = 0; i < end_effector_names.size(); ++i)
            {
                end_effector_names_array[i] = end_effector_names[i];
            }

            robot_ = std::make_shared<LeggedBody>(model_file_location, end_effector_names.size(), end_effector_names_array);
            states_ = robot_->si;
            model_file_location_ = model_file_location;
        }

        void LeggedInterface::LoadParameters(std::string parameter_file_location)
        {
            // Load the parameters from the given parameter file.
            std::map<std::string, std::string> imported_vars;
            std::ifstream file(parameter_file_location);
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

            // opts_["ipopt.fixed_variable_treatment"] = imported_vars["ipopt.fixed_variable_treatment"];
            // opts_["ipopt.max_iter"] = std::stoi(imported_vars["ipopt.max_iter"]);

            opts_["ipopt.linear_solver"] = "ma97";
            opts_["ipopt.ma97_order"] = "metis";
            // opts_["snopt.Major iterations limit"] = 1;
            // opts_["snopt.Minor iterations limit"] = 1;
            // opts_["snopt.Iterations limit"] = 1;
            // opts_["snopt.Cold Start/Warm Start"] = "Warm";

            // Extract the string value from imported_vars
            Eigen::VectorXd Q_diag = ReadVector(imported_vars["Q_diag"]);
            Eigen::VectorXd R_diag = ReadVector(imported_vars["R_diag"]);

            assert(Q_diag.size() == states_->ndx);
            assert(R_diag.size() == states_->nu);

            cost_params_.Q_diag = Q_diag;
            cost_params_.R_diag = R_diag;

            parameters_set_ = true;
        }

        void LeggedInterface::CreateProblemData(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            assert(robot_ != nullptr); // Model must be loaded
            casadi::Function L, Phi;

            CreateCost(initial_state, target_state, L, Phi);

            std::shared_ptr<opt::GeneralProblemData> gp_data = std::make_shared<opt::GeneralProblemData>(robot_->fint, robot_->fdiff, L, Phi);

            problem_data_ = std::make_shared<LeggedRobotProblemData>(gp_data,
                                                                     surfaces_,
                                                                     robot_->contact_sequence,
                                                                     states_, std::make_shared<legged::ADModel>(robot_->cmodel),
                                                                     std::make_shared<legged::ADData>(robot_->cdata),
                                                                     robot_->getEndEffectors(),
                                                                     robot_->cx, robot_->cu, robot_->cdt, initial_state, target_state);
        }

        // TODO: Generate a reference trajectory somewhere and share it betwen the objective and initial guess
        void LeggedInterface::CreateCost(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state, casadi::Function &L, casadi::Function &Phi)
        {
            assert(robot_ != nullptr);

            casadi::DM X0 = initial_state;

            assert(cost_params_.Q_diag.size() == states_->ndx);
            Eigen::MatrixXd Q_mat = cost_params_.Q_diag.asDiagonal();

            assert(cost_params_.R_diag.size() == states_->nu);
            Eigen::MatrixXd R_mat = cost_params_.R_diag.asDiagonal();

            casadi::SX Q = casadi::SX::zeros(states_->ndx, states_->ndx);
            casadi::SX R = casadi::SX::zeros(states_->nu, states_->nu);

            pinocchio::casadi::copy(Q_mat, Q);
            pinocchio::casadi::copy(R_mat, R);

            /*Both f_state_error and fdiff work for the quaternion error, but fdiff uses the quaternion logarithm, so it is more accurate albeit slightly more expensive*/
            casadi::SX X_error = robot_->f_state_error(casadi::SXVector{robot_->cx, target_state}).at(0);
            // casadi::SX X_error = robot_->fdiff(casadi::SXVector{rrobot_->cx, target_state, 1.}).at(0);

            pinocchio::computeTotalMass(robot_->cmodel, robot_->cdata);

            casadi::SX U_ref = casadi::SX::zeros(states_->nu, 1);
            // Hardcoded for 3-joint quadrupeds for now
            U_ref(casadi::Slice(2, 11, 3)) = 9.81 * robot_->cdata.mass[0] / robot_->num_end_effectors_;
            casadi::SX u_error = robot_->cu - U_ref;

            L = casadi::Function("L",
                                 {robot_->cx, robot_->cu},
                                 {1. * 0.5 * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error)) +
                                  1. * 0.5 * casadi::SX::dot(u_error, casadi::SX::mtimes(R, u_error))});

            Phi = casadi::Function("Phi",
                                   {robot_->cx},
                                   {1e3 * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error))});
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
            trajectory_opt_ = std::make_shared<LeggedTrajOpt>(problem_data_, robot_->contact_sequence, constraint_builders, decision_builder_, opts_, "ipopt");
        }

        void LeggedInterface::Update(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            UpdateProblemBoundaries(initial_state, target_state);

            // Solve the problem
            std::lock_guard<std::mutex> lock(trajectory_opt_mutex_);
            trajectory_opt_->optimize();

            std::lock_guard<std::mutex> lock(solution_mutex_);
            //@todo Akshay5312, reevaluate thread safety
            solution_interface_->UpdateSolution(trajectory_opt_->getSolutionSegments());
            solution_interface_->UpdateConstraints(trajectory_opt_->getConstraintDataSegments());
        }

        void LeggedInterface::GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const
        {
            solution_interface_->GetSolution(query_times, state_result, input_result);
        }

        void LeggedInterface::VisualizeSolutionAndConstraints(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const
        {
            std::vector<std::vector<galileo::opt::constraint_evaluations_t>> constraints = solution_interface_->GetConstraints(query_times, state_result, input_result);
            std::cout << "Size of constraints: " << constraints.size() << std::endl;

            Eigen::MatrixXd subMatrix = state_result.block(states_->q_index, 0, states_->nq, state_result.cols());
            meshcat_interface->WriteTimes(query_times, "sol_times.csv");
            meshcat_interface->WriteJointPositions(subMatrix, "sol_states.csv");
            meshcat_interface->WriteMetadata(model_file_location_, "metadata.csv");

            opt::solution::solution_t solution(query_times, state_result, input_result);
            // plotting_interface->PlotSolution(solution, {},{}, {}, {}, {}, {});
            std::cout << "Plotting constraints" << std::endl;
            plotting_interface->PlotConstraints(constraints);
        }

        void LeggedInterface::UpdateProblemBoundaries(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            std::lock_guard<std::mutex> lock(trajectory_opt_mutex_);
            problem_data_->legged_decision_problem_data.X0 = (initial_state);

            // Create a new cost
            casadi::Function L, Phi;
            CreateCost(initial_state, target_state, L, Phi);

            problem_data_->gp_data->L = L;
            problem_data_->gp_data->Phi = Phi;

            trajectory_opt_->initFiniteElements(1, initial_state);
        }

    }
}
