#include "galileo/opt/Constraint.h"
#include "galileo/legged-model/ContactSequence.h"
#include "galileo/legged-model/LeggedRobotStates.h"

namespace galileo
{
    namespace legged
    {
        namespace constraints
        {
            struct ContactConstraintProblemData
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
            };

            template <class ProblemData>
            class ContactConstraintBuilder : public opt::ConstraintBuilder<ProblemData>
            {

            public:
                ContactConstraintBuilder() : opt::ConstraintBuilder<ProblemData>() {}

                void BuildConstraint(const ProblemData &problem_data, int phase_index, opt::ConstraintData &constraint_data)
                {
                    auto mode = getModeAtPhase(problem_data, phase_index);

                    casadi::SXVector G_vec;
                    casadi::SXVector upper_bound_vec;
                    casadi::SXVector lower_bound_vec;
                    casadi::SXVector sx_foot_placement;

                    for (auto ee : problem_data.contact_constraint_problem_data.robot_end_effectors)
                    {
                        if (mode[(*ee.second)])
                        {
                            environment::SurfaceID surface = mode.getSurfaceID((*ee.second));

                            auto surface_data = (*problem_data.contact_constraint_problem_data.environment_surfaces)[surface];

                            Eigen::MatrixXd A = surface_data.A;
                            Eigen::VectorXd b = surface_data.b;
                            Eigen::VectorXd lower_bound_region_constraint = Eigen::VectorXd::Constant(b.rows(), -std::numeric_limits<double>::infinity());

                            double height = surface_data.origin_z_offset;

                            //@todo (akshay) : change this to data on num_DOF instead.

                            // auto foot_pos = casadi::SX::sym("foot_pos", 6);
                            pinocchio::SE3Tpl<galileo::opt::ADScalar, 0> frame_omf_data = problem_data.contact_constraint_problem_data.ad_data->oMf[ee.second->frame_id];

                            auto foot_pos = frame_omf_data.translation();

                            casadi::SX cfootpos = casadi::SX::sym("foot_pos", 3);
                            pinocchio::casadi::copy(foot_pos, cfootpos);

                            casadi::SX symbolic_A = casadi::SX(casadi::Sparsity::dense(A.rows(), 1));

                            pinocchio::casadi::copy(A, symbolic_A);
                            casadi::SX evaluated_vector = casadi::SX::mtimes(symbolic_A, cfootpos(casadi::Slice(1, 2)));

                            casadi::SX symbolic_b = casadi::SX(b.rows(), 1);
                            casadi::SX symbolic_lower_bound_region_constraint = casadi::SX(lower_bound_region_constraint.rows(), 1);
                            pinocchio::casadi::copy(b, symbolic_b);
                            pinocchio::casadi::copy(lower_bound_region_constraint, symbolic_lower_bound_region_constraint);
                            casadi::SX symbolic_height = height;
                            casadi::SX upper_bound = casadi::SX::vertcat({symbolic_b, symbolic_height});
                            casadi::SX lower_bound = casadi::SX::vertcat({symbolic_lower_bound_region_constraint, symbolic_height});

                            G_vec.push_back(casadi::SX::vertcat({evaluated_vector, cfootpos(0)}));
                            lower_bound_vec.push_back(lower_bound);
                            upper_bound_vec.push_back(upper_bound);

                            sx_foot_placement.push_back(cfootpos);
                        }
                    }

                    constraint_data.G = casadi::Function("G_Contact",
                                                         casadi::SXVector{problem_data.contact_constraint_problem_data.x, problem_data.contact_constraint_problem_data.u}, casadi::SXVector{casadi::SX::vertcat(G_vec)});

                    constraint_data.lower_bound = casadi::Function("lower_bound_Contact",
                                                                   casadi::SXVector{problem_data.contact_constraint_problem_data.t},
                                                                   casadi::SXVector{casadi::SX::vertcat(lower_bound_vec)});

                    constraint_data.lower_bound = casadi::Function("upper_bound_Contact",
                                                                   casadi::SXVector{problem_data.contact_constraint_problem_data.t},
                                                                   casadi::SXVector{casadi::SX::vertcat(upper_bound_vec)});

                }

            private:

                void CreateFunction(const ProblemData &problem_data, int phase_index, casadi::Function &G) const
                {
                }

                void CreateBounds(const ProblemData &problem_data, int phase_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const
                {
                }

                /**
                 * @brief getModeAtKnot gets the contact mode at the current phase
                 */
                const contact::ContactMode getModeAtPhase(const ProblemData &problem_data, int phase_index);
            };

            template <class ProblemData>
            const contact::ContactMode ContactConstraintBuilder<ProblemData>::getModeAtPhase(const ProblemData &problem_data, int phase_index)
            {
                assert(problem_data.contact_constraint_problem_data.contact_sequence != nullptr);
                return problem_data.contact_constraint_problem_data.contact_sequence->getPhase(phase_index).mode;
            }
        }
    }
}