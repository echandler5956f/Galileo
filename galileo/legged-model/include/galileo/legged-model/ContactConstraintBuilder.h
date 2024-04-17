#pragma once

#include "galileo/opt/ProblemData.h"
#include "galileo/legged-model/ContactSequence.h"
#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/centroidal.hpp>

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
                std::shared_ptr<legged::LeggedRobotStates> states;
                std::shared_ptr<legged::ADModel> ad_model;
                std::shared_ptr<legged::ADData> ad_data;
                contact::RobotEndEffectors robot_end_effectors;
                casadi::SX x;
                casadi::SX u;
                casadi::SX t;
            };

            template <class ProblemData>
            class ContactConstraintBuilder : public opt::ConstraintBuilder<ProblemData>
            {

            public:
                ContactConstraintBuilder() : opt::ConstraintBuilder<ProblemData>() {}

                void buildConstraint(const ProblemData &problem_data, int phase_index, opt::ConstraintData &constraint_data)
                {
                    auto mode = getModeAtPhase(problem_data, phase_index);

                    casadi::SXVector G_vec;
                    casadi::SXVector upper_bound_vec;
                    casadi::SXVector lower_bound_vec;
                    for (auto ee : problem_data.contact_constraint_problem_data.robot_end_effectors)
                    {
                        if (mode[(*ee.second)])
                        {
                            environment::SurfaceID surface = mode.getSurfaceID((*ee.second));

                            environment::SurfaceData surface_data = (*problem_data.contact_constraint_problem_data.environment_surfaces)[surface];

                            // Get surface data, Copy into symbolic Casadi Expr.
                            int n_constraints = surface_data.A.rows();
                            casadi::SX symbolic_A = casadi::SX(n_constraints, 2);
                            pinocchio::casadi::copy(surface_data.A, symbolic_A);

                            casadi::SX symbolic_b = casadi::SX(n_constraints, 1);
                            pinocchio::casadi::copy(surface_data.b, symbolic_b);

                            casadi::SX symbolic_surface_translation(3, 1);
                            pinocchio::casadi::copy(surface_data.surface_transform.translation(), symbolic_surface_translation);

                            casadi::SX symbolic_surface_rotation(3, 3);
                            pinocchio::casadi::copy(surface_data.surface_transform.rotation(), symbolic_surface_rotation);

                            // Defined Bounds
                            casadi::SX lower_bound = casadi::SX(n_constraints, 1);
                            pinocchio::casadi::copy(
                                Eigen::VectorXd::Constant(n_constraints, -std::numeric_limits<double>::infinity()),
                                lower_bound);

                            casadi::SX upper_bound = casadi::SX(n_constraints, 1);
                            pinocchio::casadi::copy(
                                Eigen::VectorXd::Constant(n_constraints, 0),
                                upper_bound);

                            // Get foot position in global frame
                            pinocchio::SE3Tpl<legged::ADScalar, 0> frame_omf_data = problem_data.contact_constraint_problem_data.ad_data->oMf[ee.first];
                            auto foot_pos = frame_omf_data.translation();
                            casadi::SX c_foot_pos_in_world = casadi::SX::sym("foot_pos", 3, 1);
                            pinocchio::casadi::copy(foot_pos, c_foot_pos_in_world);

                            // Get foot position in surface frame
                            casadi::SX foot_pos_offset = (c_foot_pos_in_world - symbolic_surface_translation);
                            casadi::SX foot_pos_in_surface = casadi::SX::mtimes(symbolic_surface_rotation, foot_pos_offset);

                            casadi::SX evaluated_vector = casadi::SX::mtimes(symbolic_A, foot_pos_in_surface(casadi::Slice(0, 2))) - symbolic_b;

                            // // if all zeros, then infinite ground
                            // if (surface_data.A(0) != 0 && surface_data.A(1) != 0)
                            // {
                            //     std::cout << "HEREREWFWFW" << std::endl;
                            //     G_vec.push_back(evaluated_vector);

                            //     lower_bound_vec.push_back(lower_bound);
                            //     upper_bound_vec.push_back(upper_bound);
                            // }

                            G_vec.push_back(foot_pos_in_surface(2));
                            lower_bound_vec.push_back(casadi::SX::zeros(1, 1));
                            upper_bound_vec.push_back(casadi::SX::zeros(1, 1));
                        }
                    }

                    constraint_data.G = casadi::Function("G_Contact",
                                                         casadi::SXVector{problem_data.contact_constraint_problem_data.x, problem_data.contact_constraint_problem_data.u}, casadi::SXVector{vertcat(G_vec)});

                    constraint_data.lower_bound = casadi::Function("lower_bound_Contact",
                                                                   casadi::SXVector{problem_data.contact_constraint_problem_data.t},
                                                                   casadi::SXVector{vertcat(lower_bound_vec)});

                    constraint_data.upper_bound = casadi::Function("upper_bound_Contact",
                                                                   casadi::SXVector{problem_data.contact_constraint_problem_data.t},
                                                                   casadi::SXVector{vertcat(upper_bound_vec)});

                    constraint_data.metadata.name = "Contact Constraint";
                }

            private:
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