#include "variables/include/Constraint.h"

namespace acro
{
    namespace model
    {
        namespace legged
        {

            struct FrictionConeProblemData
            {

                std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;
                std::shared_ptr<contact::ContactSequence> contact_sequence;

                std::shared_ptr<pinocchio::Model> model;
                RobotEndEffectors robot_end_effectors;

                int num_knots;
                int collocation_points_per_knot;

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
             *      or      [0   0  mu] * Rotation * GRF
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
            class FrictionConeConstraintBuilder : Constraint<ProblemData>
            {

            public:
                FrictionConeConstraintBuilder() : ConstraintBuilder() {}

                /**
                 *
                 * @brief Generate flags for each knot point. We set it to all ones, applicable at each knot.
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "FrictionConeProblemData" NAMED "friction_cone_problem_data"
                 * @param apply_at
                 */
                void CreateApplyAt(const ProblemData &problem_data, int knot_index, Eigen::VectorXi &apply_at) const override
                {
                    uint num_points = problem_data.friction_cone_problem_data.collocation_points_per_knot;
                    apply_at = Eigen::VectorXi::Constant(num_points, 1);
                }

                /**
                 * @brief Generate bounds for a vector of points.
                 *
                 * For both approximations, each value is less than 0, with no lower bound.
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "FrictionConeProblemData" NAMED "friction_cone_problem_data"
                 * @param upper_bound
                 * @param lower_bound
                 */
                void CreateBounds(const ProblemData &problem_data, int knot_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const override;

                /**
                 * @brief Generate a function to evaluate each point
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "FrictionConeProblemData" NAMED "friction_cone_problem_data"
                 * @param G
                 */
                void CreateFunction(const ProblemData &problem_data, int knot_index, casadi::Function &G) const;

            private:
                uint getNumConstraintPerEEPerState(const ProblemData &problem_data) const;

                void CreateSingleEndEffectorFunction(const std::string &EndEffectorID, const ProblemData &problem_data, int knot_index, casadi::Function &G) const;

                const contact::ContactMode &getModeAtKnot(const ProblemData &problem_data, int knot_index);

                Eigen::Matrix<double, 3, 3> getContactSurfaceRotationAtMode(const std::string &EndEffectorID, const ProblemData &problem_data, const contact::ContactMode &mode);

                Eigen::MatrixXd getConeConstraintApproximation(const ProblemData &problem_data);
            };

            template <class ProblemData>
            void FrictionConeConstraintBuilder<ProblemData>::CreateBounds(const ProblemData &problem_data, casadi::Function &upper_bound, casadi::Function &lower_bound) const override;
            {
                uint num_points = problem_data.friction_cone_problem_data.collocation_points_per_knot;

                uint num_constraints = num_points * getNumConstraintPerEEPerState(problem_data);

                Eigen::VectorXd upper_bound_vect = Eigen::VectorXd::Constant(num_constraints, 0);
                Eigen::VectorXd lower_bound_vect = Eigen::VectorXd::Constant(num_constraints, -std::numeric_limits<double>::infinity());

                // @TODO create a casadi::Function to return vectors
            }

            template <class ProblemData>
            void FrictionConeConstraintBuilder<ProblemData>::CreateFunction(const ProblemData &problem_data, casadi::Function &G) const{
                //CREATE CASADI MAP FROM EACH END EFFECTOR 
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
            void FrictionConeConstraintBuilder<ProblemData>::CreateSingleEndEffectorFunction(const std::string &EndEffectorID, const ProblemData &problem_data, casadi::Function &G) const
            {
                // Get the size of the constraint applied to an end effector at a state.
                uint num_contraint_per_ee_per_point = getNumConstraintPerEEPerState(problem_data);

                // Get the mode at the knot point.
                const contact::ContactMode &mode = getModeAtKnot(problem_data, knot_index);

                // If the end effector is not in contact, Do not apply the constraint; the function is all zeros
                if (!mode.combination_definition[EndEffectorID])
                {
                    // SET ALL ZEROS CASADI FUNCTION
                    //  Function = Eigen::MatrixXd::Zero(num_contraint_per_ee_per_point, 3) * GRF
                    return;
                }

                Eigen::Matrix<double, 3, 3> rotation = getContactSurfaceRotationAtMode(EndEffectorID, problem_data, mode);
            
                Eigen::MatrixXd cone_constraint_approximation = getConeConstraintApproximation(problem_data);

                Eigen::MatrixXd rotated_cone_constraint = cone_constraint_approximation * rotation;

                FrictionConeProblemData::ApproximationOrder approximation_order = problem_data.friction_cone_problem_data.approximation_order;

                if (approximation_order == FrictionConeProblemData::ApproximationOrder::FIRST_ORDER)
                {
                    //SET CASADI FUNCTION 
                    // Function = rotated_cone_constraint * GRF(EndEffector)
                }
                else
                {
                    //SET CASADI FUNCTION 
                    // evaluated_vector = (rotated_cone_constraint * GRF(EndEffector))
                    // Function = (evaluated_vector[0] - evaluated_vector.tail(2).normSquared());
                }
            }

            template <class ProblemData>
            const contact::ContactMode &FrictionConeConstraintBuilder<ProblemData>::getModeAtKnot(const ProblemData &problem_data, int knot_index)
            {
                assert(problem_data.contact_sequence != nullptr);

                contact::ContactSequence::Phase phase;
                contact::ContactSequence::CONTACT_SEQUENCE_ERROR status;
                problem_data.contact_sequence->getPhaseAtKnot(knot_index, phase, status);

                assert(status == contact::ContactSequence::CONTACT_SEQUENCE_ERROR::OK);

                return phase.mode;
            }

            template <class ProblemData>
            Eigen::Matrix<double, 3, 3> FrictionConeConstraintBuilder<ProblemData>::getContactSurfaceRotationAtMode(const std::string &EndEffectorID, const ProblemData &problem_data, const contact::ContactMode &mode)
            {
                assert(problem_data.environment_surfaces != nullptr);

                contact::SurfaceID contact_surface_id = mode.getSurfaceIDForEE(EndEffectorID);

                if (contact_surface_id == contact::NO_SURFACE)
                {
                    // This is undesired behavior; for now we just return a matrix of zeros.
                    return Eigen::Matrix<double, 3, 3>::Zero();
                }

                return (*problem_data.environment_surfaces)[contact_surface_id].Rotation();
            }


            template <class ProblemData>
            Eigen::MatrixXd FrictionConeConstraintBuilder<ProblemData>::getConeConstraintApproximation(const ProblemData &problem_data)
            {
                Eigen::MatrixXd cone_constraint_approximation(num_contraint_per_ee_per_point, 3);
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
                    cone_constraint_approximation.block(0, 0, 2, 2) = Eigen::MatrixXd::Identity(2);
                    cone_constraint_approximation.block(2, 0, 2, 2) = -Eigen::MatrixXd::Identity(2);

                    cone_constraint_approximation.block(0, 2, 4, 1) = Eigen::MatrixXd::VectorXd::Constant(4, -mu);

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
                    cone_constraint_approximation.block(1, 0, 2, 2) = Eigen::MatrixXd::Identity(2);
                }

                return cone_constraint_approximation;
            }

        }
    }
}