#include "galileo/legged-model/LeggedInterface.h"

namespace galileo
{
    namespace legged
    {
        LeggedInterface::LeggedInterface()
        {
            surfaces_ = std::make_shared<galileo::legged::environment::EnvironmentSurfaces>();
            solution_interface_ = std::make_shared<galileo::opt::solution::Solution>();
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
        }

        void LeggedInterface::LoadParameters(std::string parameter_file_location)
        {
            // Hardcoded parameters for now
            // Load the parameters from the given parameter file.
            opts_["ipopt.linear_solver"] = "ma97";
            opts_["ipopt.ma97_order"] = "metis";
            opts_["ipopt.fixed_variable_treatment"] = "make_constraint";
            opts_["ipopt.max_iter"] = 250;

            cost_params_.R_diag = Eigen::VectorXd(states_->nu);
            cost_params_.R_diag << 1e-3, 1e-3, 1e-3,                        /*First contact wrench error weights*/
                1e-3, 1e-3, 1e-3,                                           /*Second contact wrench error weights*/
                1e-3, 1e-3, 1e-3,                                           /*Third contact wrench error weights*/
                1e-3, 1e-3, 1e-3,                                           /*Fourth contact wrench error weights*/
                10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.; /*Joint acceleration error weights*/

            cost_params_.Q_diag = Eigen::VectorXd(states_->ndx);
            cost_params_.Q_diag << 15., 15., 30., 5., 10., 10.,             /*Centroidal momentum error weights*/
                0., 0., 0., 0., 0., 0.,                                     /*Rate of Centroidal momentum error weights*/
                500., 500., 500., 0.1, 0.1, 0.1,                            /*Floating base position and orientation (exponential coordinates) error weights*/
                20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., /*Joint position error weights*/
                0., 0., 0., 0., 0., 0.,                                     /*Floating base velocity error weights*/
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;             /*Joint velocity error weights*/
        }

        void LeggedInterface::CreateProblemData(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
            casadi::Function L, Phi;

            CreateCost(initial_state, target_state, L, Phi);

            std::shared_ptr<opt::GeneralProblemData> gp_data = std::make_shared<opt::GeneralProblemData>(robot_->fint, robot_->fdiff, L, Phi);

            problem_data_ = std::make_shared<LeggedRobotProblemData>(gp_data,
                                                                     surfaces_,
                                                                     robot_->contact_sequence,
                                                                     states_, std::make_shared<galileo::opt::ADModel>(robot_->cmodel),
                                                                     std::make_shared<galileo::opt::ADData>(robot_->cdata),
                                                                     robot_->getEndEffectors(),
                                                                     robot_->cx, robot_->cu, robot_->cdt, initial_state);
        }

        // TODO: Generate a reference trajectory somewhere and share it betwen the objective and initial guess
        void LeggedInterface::CreateCost(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state, casadi::Function &L, casadi::Function &Phi)
        {
            // Hardcoded costs for now
            // Updates the robot running and terminal costs L and Phi.

            // HARDCODE COST WEIGHTS
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
                                   {1. * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error))});
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
            // Hardcoded initialization for now
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

            trajectory_opt_ = std::make_shared<LeggedTrajOpt>(problem_data_, robot_->contact_sequence, constraint_builders, decision_builder_, opts_);
        }

        void LeggedInterface::Update(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state, casadi::MXVector &solution, Eigen::VectorXd &times)
        {
            UpdateProblemBoundaries(initial_state, target_state);

            // Solve the problem
            solution = trajectory_opt_->optimize();
            solution_segments_.clear();
            trajectory_opt_->getSolutionSegments(solution_segments_, times);
            solution_interface_->UpdateSolution(solution_segments_);
            solution_ = solution;
        }

        void LeggedInterface::UpdateProblemBoundaries(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state)
        {
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
