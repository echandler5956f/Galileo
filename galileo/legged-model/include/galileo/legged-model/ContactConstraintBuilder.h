#include "galileo/variables/Constraint.h"
#include "galileo/model/ContactSequence.h"
#include "galileo/variables/States.h"

namespace galileo
{
    namespace model
    {
        namespace legged
        {

            struct ContactConstraintProblemData
            {

                std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;
                std::shared_ptr<contact::ContactSequence> contact_sequence;
                std::shared_ptr<variables::States> states;
                std::shared_ptr<pinocchio::Model> model;
                RobotEndEffectors robot_end_effectors;
                casadi::SX x; // this needs to be initialized to casadi::SX::sym("x", states->nx) somewhere
                casadi::SX u; // this needs to be initialized to casadi::SX::sym("u", states->nu) somewhere
                casadi::SX t; // this needs to be initialized to casadi::SX::sym("t") somewhere
                int num_knots;
                int collocation_points_per_knot;
            };

            template <class ProblemData>
            class ContactConstraintProblemData : Constraint<ProblemData>
            {

            public:
                ContactConstraintProblemData() : ConstraintBuilder() {}

                void BuildConstraint(const ProblemData &problem_data, int knot_index, ConstraintData &constraint_data)
                {
                    auto mode = getModeAtKnot(problem_data, knot_index);

                    for (auto ee : problem_data.robot_end_effectors)
                    {
                        if (mode[(*ee.second)])
                        {
                            environment::SurfaceID surface = mode.getSurfaceID((*ee.second));

                            auto surface_data = (*problem_data.environment_surfaces)[surface];

                            Eigen::MatrixXd A = surface_data.A;
                            Eigen::VectorXd b = surface_data.b;

                            double height = surface_data.height;

                            // @todo(ethan) : add the constraints (height = ee_position.bottom(1)) & (A * ee_position.top(2) - b) <= 0 to the constraint data

                            // Function FX(x, end_effector)
                        }
                    }
                }

            private:
                /**
                 *
                 * @brief Generate flags for each knot point. We set it to all ones, applicable at each knot.
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "ContactConstraintProblemData" NAMED "contact_constraint_problem_data"
                 * @param apply_at
                 */
                void CreateApplyAt(const ProblemData &problem_data, int knot_index, Eigen::VectorXi &apply_at) const override
                {
                    uint num_points = problem_data.contact_constraint_problem_data.collocation_points_per_knot;
                    apply_at = Eigen::VectorXi::Constant(num_points, 1);
                }

                /**
                 * @brief getModeAtKnot gets the contact mode at the current knot
                 */
                const contact::ContactMode &getModeAtKnot(const ProblemData &problem_data, int knot_index);
            };

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
        }
    }
}