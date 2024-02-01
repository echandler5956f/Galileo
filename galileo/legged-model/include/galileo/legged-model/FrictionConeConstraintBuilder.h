#include "galileo/opt/Constraint.h"
#include "galileo/legged-model/ContactSequence.h"
#include "galileo/legged-model/LeggedRobotStates.h"

namespace galileo
{
    namespace legged
    {
        namespace constraints
        {

            struct FrictionConeProblemData
            {
                std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;
                std::shared_ptr<contact::ContactSequence> contact_sequence;
                std::shared_ptr<opt::States> states;
                std::shared_ptr<opt::ADModel> ad_model;
                std::shared_ptr<opt::ADData> ad_data;
                contact::RobotEndEffectors robot_end_effectors;
                casadi::SX x;
                casadi::SX u;
                casadi::SX t;
                int num_knots;

                float mu;

                enum ApproximationOrder
                {
                    FIRST_ORDER,
                    SECOND_ORDER
                };
                ApproximationOrder approximation_order;
            };

            /**
             * A builder for the Friction Cone Constraint.
             *
             * First Order Approximation --
             *      [ 1    0  -mu]
             *      [ 0    1  -mu]
             *      [-1    0  -mu] * Rotation * GRF  <= 0
             *      [ 0   -1  -mu]
             *      [ 0    0   -1]
             *
             *
             * Second Order Approximation --
             *                              [0   0  mu]
             *      LorentzConeConstraint ( [0   1   0] * Rotation * GRF )
             *                              [1   0   0]
             *
             *      or     [0   0  mu] * Rotation * GRF
             *               - 2norm(   [0   1   0] * Rotation * GRF  )  => 0
             *                          [1   0   0]
             *
             * "Rotation" must be a transformation such that (Rotation * unit_surface_normal = [0 0 1]).
             *
             * Applied at each knot point.
             *
             * @tparam ProblemData must contain an instance of "FrictionConeProblemData" named "friction_cone_problem_data"
             */
            template <class ProblemData>
            class FrictionConeConstraintBuilder : public opt::ConstraintBuilder<ProblemData>
            {

            public:
                FrictionConeConstraintBuilder() : opt::ConstraintBuilder<ProblemData>() {}

            private:
                /**
                 * @brief Generate bounds for a vector of points.
                 *
                 * For both approximations, each value is less than 0, with no lower bound.
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "FrictionConeProblemData" NAMED "friction_cone_problem_data"
                 * @param phase_index the index of the hase
                 * @param upper_bound Lower bound of the constraint at each point
                 * @param lower_bound Upper bound of the constraint at each point
                 */
                void CreateBounds(const ProblemData &problem_data, int phase_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const;

                /**
                 * @brief Generate a function to evaluate each point
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "FrictionConeProblemData" NAMED "friction_cone_problem_data"
                 * @param phase_index the index of the phase
                 * @param G A function that evaluates the constraint at each point
                 */
                void CreateFunction(const ProblemData &problem_data, int phase_index, casadi::Function &G) const;

                /**
                 * @brief get the number of constraints applied to each end effector at a given state. 5 for a first order approximation, 1 for second order
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "FrictionConeProblemData" NAMED "friction_cone_problem_data"
                 * @return uint
                 */
                uint getNumConstraintPerEEPerState(const ProblemData &problem_data) const;

                /**
                 * @brief Get the Casadi::Function object, G, to be applied to the end effector given by "EndEffectorID" at a knot point.
                 * For a Second Order Constraint, for instance, this is a function that returns 11 value, the result of the lorentz cone constraint.
                 * function = g( state ) @ EndEffectorID
                 */
                void CreateSingleEndEffectorFunction(const std::string &EndEffectorID, const ProblemData &problem_data, int phase_index, const casadi::SX &u_in, casadi::SX &G_out) const;

                /**
                 * @brief getModeAtKnot gets the contact mode at the current phase
                 */
                const contact::ContactMode getModeAtPhase(const ProblemData &problem_data, int phase_index) const;

                /**
                 * @brief getContactSurfaceRotationAtMode gets a rotation matrix describing the normal to the contact surface at the current knot.
                 * If the End Effector is in contact with a surface, the result is a rotation describing the surface noraml, else it is a matrix of zeros.
                 */
                Eigen::Matrix<double, 3, 3> getContactSurfaceRotationAtMode(const std::string &EndEffectorID, const ProblemData &problem_data, const contact::ContactMode &mode) const;

                /**
                 * @brief Each approximated cone constraint can be represented as an operation on a transformed "ground reaction force".
                 *           For a First order approximation, the constraint is of the form A_first_order * rotated_ground_reation_force <= 0.
                 *           For a second order approximation, this is of the form LorentzConeConstraint(A_second_order * rotated_ground_reation_force) <= 0
                 *           This function returns "A".
                 */
                Eigen::MatrixXd getConeConstraintApproximation(const ProblemData &problem_data) const;
            };

            template <class ProblemData>
            void FrictionConeConstraintBuilder<ProblemData>::CreateBounds(const ProblemData &problem_data, int phase_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const
            {
                uint num_points = problem_data.friction_cone_problem_data.num_knots;

                uint num_constraints = num_points * getNumConstraintPerEEPerState(problem_data);

                Eigen::VectorXd upper_bound_vect = Eigen::VectorXd::Constant(num_constraints, 0);
                Eigen::VectorXd lower_bound_vect = Eigen::VectorXd::Constant(num_constraints, -std::numeric_limits<double>::infinity());

                // Convert Eigen to casadi (see traj_test example)
                casadi::SX upper_bound_casadi = casadi::SX(casadi::Sparsity::dense(upper_bound_vect.rows(), 1));
                pinocchio::casadi::copy(upper_bound_vect, upper_bound_casadi);

                casadi::SX lower_bound_casadi = casadi::SX(casadi::Sparsity::dense(lower_bound_vect.rows(), 1));
                pinocchio::casadi::copy(lower_bound_vect, lower_bound_casadi);

                upper_bound = casadi::Function("upper_bound", casadi::SXVector{problem_data.friction_cone_problem_data.t}, casadi::SXVector{upper_bound_casadi});
                lower_bound = casadi::Function("lower_bound", casadi::SXVector{problem_data.friction_cone_problem_data.t}, casadi::SXVector{lower_bound_casadi});
            }

            template <class ProblemData>
            void FrictionConeConstraintBuilder<ProblemData>::CreateFunction(const ProblemData &problem_data, int phase_index, casadi::Function &G) const
            {
                // CREATE CASADI MAP FROM EACH END EFFECTOR
                // CreateSingleEndEffectorFunction creates a function, g(state) @ end_effector. Here, we must get the constraint for each end_effector at this knot point, and apply them to each collocation point in the knot.
                // Possibly, this would mean creating a map.
                // In doing so, we create a a function that evaluates each end effector at each collocation point in the knot segment.
                casadi::SXVector G_vec;
                casadi::SX u_in = casadi::SX::sym( "u", problem_data.friction_cone_problem_data.states->nu() );
                for (auto &end_effector : problem_data.friction_cone_problem_data.robot_end_effectors)
                {
                    casadi::SX G_out;
                    CreateSingleEndEffectorFunction(end_effector.first, problem_data, phase_index, u_in, G_out);
                    if (!G_out.is_zero())
                    {
                        G_vec.push_back(G_out);
                        u_vec.push_back(u_in);
                    }
                }
                G = casadi::Function("G_FrictionCone", casadi::SXVector{problem_data.friction_cone_problem_data.x, casadi::SX::vertcat(u_vec)}, casadi::SXVector{casadi::SX::vertcat(G_vec)});
            }

            template <class ProblemData>
            uint FrictionConeConstraintBuilder<ProblemData>::getNumConstraintPerEEPerState(const ProblemData &problem_data) const
            {
                FrictionConeProblemData::ApproximationOrder approximation_order = problem_data.friction_cone_problem_data.approximation_order;
                if (approximation_order == FrictionConeProblemData::ApproximationOrder::SECOND_ORDER)
                {
                    return 1;
                }
                return 5;
            }

            template <class ProblemData>
            void FrictionConeConstraintBuilder<ProblemData>::CreateSingleEndEffectorFunction(const std::string &EndEffectorID, const ProblemData &problem_data, int phase_index, const casadi::SX &u_in, casadi::SX &G_out) const
            {
                // Get the size of the constraint applied to an end effector at a state.
                uint num_contraint_per_ee_per_point = getNumConstraintPerEEPerState(problem_data);

                // Get the mode at the knot point.
                const contact::ContactMode mode = getModeAtPhase(problem_data, phase_index);
                casadi::SX x = problem_data.friction_cone_problem_data.x;
                auto it = mode.combination_definition.find(EndEffectorID);
                bool dof6 = (it != mode.combination_definition.end()) ? problem_data.friction_cone_problem_data.robot_end_effectors.find(EndEffectorID)->second->is_6d : false;

                // If the end effector is not in contact, Do not apply the constraint; the function is all zeros
                if (it != mode.combination_definition.end() && !it->second)

                {
                    // SET ALL ZEROS CASADI FUNCTION
                    //  Function = Eigen::MatrixXd::Zero(num_contraint_per_ee_per_point, 3) * GRF
                    G_out = casadi::SX(0);
                    return;
                }

                auto u_ee = problem_data.friction_cone_problem_data.states->get_f(u_in, EndEffectorID);

                Eigen::Matrix<double, 3, 3> rotation = getContactSurfaceRotationAtMode(EndEffectorID, problem_data, mode);

                Eigen::MatrixXd cone_constraint_approximation = getConeConstraintApproximation(problem_data);

                Eigen::MatrixXd rotated_cone_constraint = cone_constraint_approximation * rotation;

                FrictionConeProblemData::ApproximationOrder approximation_order = problem_data.friction_cone_problem_data.approximation_order;

                casadi::SX symbolic_rotated_cone_constraint = casadi::SX(casadi::Sparsity::dense(rotated_cone_constraint.rows(), 1));
                pinocchio::casadi::copy(rotated_cone_constraint, symbolic_rotated_cone_constraint);

                // great care must be taken in doing this to ensure that the appropriate u is passed into this function.
                auto evaluated_vector = casadi::SX::mtimes(symbolic_rotated_cone_constraint, u_ee);
                if (approximation_order == FrictionConeProblemData::ApproximationOrder::FIRST_ORDER)
                {
                    //@todo (ethan)
                    // SET CASADI FUNCTION
                    //  Function = rotated_cone_constraint * GRF(EndEffector)
                    G_out = evaluated_vector;
                }
                else
                {
                    //@todo (ethan)
                    // SET CASADI FUNCTION
                    //  evaluated_vector = (rotated_cone_constraint * GRF(EndEffector))
                    //  Function = (evaluated_vector[0] - evaluated_vector.tail(2).normSquared());
                    G_out = evaluated_vector(0) - casadi::SX::norm_2(evaluated_vector(casadi::Slice(1, 3)));
                }
            }

            template <class ProblemData>
            const contact::ContactMode FrictionConeConstraintBuilder<ProblemData>::getModeAtPhase(const ProblemData &problem_data, int phase_index) const
            {
                assert(problem_data.friction_cone_problem_data.contact_sequence != nullptr);
                return problem_data.friction_cone_problem_data.contact_sequence->getPhase(phase_index).mode;
            }

            template <class ProblemData>
            Eigen::Matrix<double, 3, 3> FrictionConeConstraintBuilder<ProblemData>::getContactSurfaceRotationAtMode(const std::string &EndEffectorID, const ProblemData &problem_data, const contact::ContactMode &mode) const
            {
                assert(problem_data.friction_cone_problem_data.environment_surfaces != nullptr);

                environment::SurfaceID contact_surface_id = mode.getSurfaceID(EndEffectorID);

                if (contact_surface_id == environment::NO_SURFACE)
                {
                    // This is undesired behavior; for now we just return a matrix of zeros.
                    return Eigen::Matrix<double, 3, 3>::Zero();
                }

                return (*problem_data.friction_cone_problem_data.environment_surfaces)[contact_surface_id].Rotation();
            }

            template <class ProblemData>
            Eigen::MatrixXd FrictionConeConstraintBuilder<ProblemData>::getConeConstraintApproximation(const ProblemData &problem_data) const
            {
                Eigen::MatrixXd cone_constraint_approximation(getNumConstraintPerEEPerState(problem_data), 3);
                FrictionConeProblemData::ApproximationOrder approximation_order = problem_data.friction_cone_problem_data.approximation_order;

                float mu = problem_data.friction_cone_problem_data.mu;
                assert(mu >= 0);

                if (approximation_order == FrictionConeProblemData::ApproximationOrder::FIRST_ORDER)
                {
                    /**
                     *  First Order Approximation --
                     *      [ 1    0  -mu]
                     *      [ 0    1  -mu]
                     *      [-1    0  -mu] * Rotation * GRF  <= 0
                     *      [ 0   -1  -mu]
                     *      [ 0    0   -1]
                     *
                     */
                    cone_constraint_approximation.block(0, 0, 2, 2) = Eigen::MatrixXd::Identity(2, 2);
                    cone_constraint_approximation.block(2, 0, 2, 2) = -Eigen::MatrixXd::Identity(2, 2);

                    cone_constraint_approximation.block(0, 2, 4, 1) = Eigen::VectorXd::Constant(4, -mu);

                    cone_constraint_approximation(4, 2) = -1;
                }
                else
                {
                    /**
                     * Second Order Approximation --
                     *                              [0   0  mu]
                     *      LorentzConeConstraint ( [0   1   0] * Rotation * GRF )
                     *                              [1   0   0]
                     */
                    cone_constraint_approximation(0, 2) = mu;
                    cone_constraint_approximation.block(1, 0, 2, 2) = Eigen::MatrixXd::Identity(2, 2);
                }

                return cone_constraint_approximation;
            }
        }
    }
}