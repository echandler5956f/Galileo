#include "galileo/legged-model/LeggedInterface.h"

namespace galileo
{
    namespace legged
    {
        LeggedInterface::LeggedInterface()
        {
            surfaces_ = std::make_shared<galileo::legged::environment::EnvironmentSurfaces>();
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
            opts_["ipopt.fixed_variable_treatment"] = "make_constraint";
            opts_["ipopt.max_iter"] = 10;

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

            std::shared_ptr<opt::GeneralProblemData> gp_data = std::make_shared<opt::GeneralProblemData>(robot_->fint, robot_->fdif, L, Phi);

            problem_data_ = std::make_shared<LeggedRobotProblemData>(gp_data,
                                                                     surfaces_,
                                                                     robot_->contact_sequence,
                                                                     states_, std::make_shared<galileo::opt::ADModel>(robot_->cmodel),
                                                                     std::make_shared<galileo::opt::ADData>(robot_->cdata),
                                                                     robot_->getEndEffectors(),
                                                                     robot_->cx, robot_->cu, robot_->cdt, initial_state);
        }

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

            auto qf = states_->get_q(target_state);
            casadi::SX target_pos = qf(casadi::Slice(0, 3));
            // Default target rot is identity
            casadi::SX target_rot = casadi::SX::eye(3);

            // GET TARGET POS AND ROT FROM THE TARGET STATE

            pinocchio::SE3Tpl<galileo::opt::ADScalar, 0> oMf = robot_->cdata.oMf[robot_->model.getFrameId(body_name_, pinocchio::BODY)];
            auto rot = oMf.rotation();
            auto pos = oMf.translation();
            casadi::SX crot = casadi::SX::zeros(3, 3);
            casadi::SX cpos = casadi::SX::zeros(3, 1);
            pinocchio::casadi::copy(rot, crot);
            pinocchio::casadi::copy(pos, cpos);

            casadi::SX rot_c = casadi::SX::mtimes(crot.T(), target_rot);

            casadi::SX target_error_casadi = casadi::SX::vertcat(casadi::SXVector{cpos - target_pos, casadi::SX::inv_skew(rot_c - rot_c.T()) / 2});

            casadi::SX X_ref = casadi::SX(X0);
            X_ref(casadi::Slice(states_->nh + states_->ndh, states_->nh + states_->ndh + 3)) = target_pos;

            auto h_and_dh_error = robot_->cx(casadi::Slice(0, states_->nh + states_->ndh)) - X_ref(casadi::Slice(0, states_->nh + states_->ndh));
            auto x_error = robot_->cx(casadi::Slice(states_->nh + states_->ndh + states_->nqb, states_->nx)) - X_ref(casadi::Slice(states_->nh + states_->ndh + states_->nqb, states_->nx));

            casadi::SX X_error = casadi::SX::vertcat({h_and_dh_error,
                                                      target_error_casadi,
                                                      x_error});

            pinocchio::computeTotalMass(robot_->model, robot_->data);

            casadi::SX U_ref = casadi::SX::zeros(states_->nu, 1);
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

            bool fully_initialized_ = true;
        }

        void LeggedInterface::CreateTrajOpt()
        {
            auto constraint_builders = getLeggedConstraintBuilders();
            decision_builder_ = std::make_shared<galileo::legged::constraints::LeggedDecisionDataBuilder<LeggedRobotProblemData>>();

            trajectory_opt_ = std::make_shared<LeggedTrajOpt>(problem_data_, robot_->contact_sequence, constraint_builders, decision_builder_, opts_);
        }

        void LeggedInterface::Update(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state, casadi::MXVector &solution)
        {
            UpdateProblemBoundaries(initial_state, target_state);

            // Solve the problem
            solution = trajectory_opt_->optimize();
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
