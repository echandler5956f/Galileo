#pragma once

#include "galileo/opt/PhaseSequence.h"
#include "galileo/legged-model/LeggedRobotStates.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"

#include <pinocchio/autodiff/casadi.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/data.hpp>
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
        namespace contact
        {

            /**
             * @brief A struct for holding the contact mode of the robot.
             *
             */
            class ContactMode
            {
            public:
                ContactMode() {}

                /**
                 * @brief Error codes
                 *
                 */
                enum ContactModeValidity
                {
                    VALID,               // No error
                    SURFACE_NOT_DEFINED, // The surface is not defined
                    DIFFERING_SIZES      // The sizes of the contact combination and contact surfaces are different
                };
                /**
                 * @brief Maps EE's to flags (true = in contact).
                 *
                 */
                ContactCombination combination_definition;

                /**
                 * @brief Gets the contact combination for a given end effector.
                 *
                 * @param ee The end effector to get the contact combination for.
                 */
                bool &operator[](const EndEffector &ee)
                {
                    return combination_definition[ee.frame_name];
                }

                /**
                 * @brief Gets the contact combination for a given end effector.
                 *
                 * @param ee_name The name of the end effector to get the contact combination for.
                 */
                bool &operator[](const std::string &ee_name)
                {
                    return combination_definition[ee_name];
                }

                const bool &at(const EndEffector &ee) const
                {
                    return combination_definition.at(ee.frame_name);
                }

                const bool &at(const std::string &ee_name) const
                {
                    return combination_definition.at(ee_name);
                }

                /**
                 * @brief Gets the contact surface for a given end effector in this mode.
                 *
                 * @param ee_name The end effector to get the contact surface for.
                 */
                const environment::SurfaceID &getSurfaceID(const std::string &ee_name) const;

                /**
                 * @brief Gets the contact surface for a given end effector in this mode.
                 *
                 * @param ee The end effector to get the contact surface for.
                 */
                const environment::SurfaceID &getSurfaceID(const EndEffector &ee) const
                {
                    return getSurfaceID(ee.frame_name);
                }

                /**
                 * @brief Gets which surfaces the EEs are in contact with.
                 *
                 */
                std::vector<environment::SurfaceID> contact_surfaces;

                /**
                 * @brief Makes the combination valid. If an EE is not in contact, it makes the corresponding contact surface NO_SURFACE.
                 *
                 * @param validity Error code
                 */
                void MakeValid(ContactModeValidity &validity);

                /**
                 * @brief creates the dynamics for this mode
                 *
                 */
                void createModeDynamics(std::shared_ptr<opt::Model> model, RobotEndEffectors end_effectors, std::shared_ptr<opt::LeggedRobotStates> states)
                {
                    opt::Data data = opt::Data(*model);

                    auto si = states;
                    casadi::SX cx = casadi::SX::sym("x", si->nx);
                    casadi::SX cu = casadi::SX::sym("u", si->nu);

                    // Get the slices of the state that represent
                    // Joint pos/vel
                    auto cq = si->get_q(cx);
                    auto cv = si->get_v(cx);

                    opt::ADModel cmodel = (*model).cast<opt::ADScalar>();
                    opt::ADData cdata(cmodel);
                    opt::ConfigVectorAD q_AD(model->nq);
                    opt::TangentVectorAD v_AD(model->nv);

                    q_AD = Eigen::Map<opt::ConfigVectorAD>(static_cast<std::vector<opt::ADScalar>>(cq).data(), model->nq, 1);
                    v_AD = Eigen::Map<opt::TangentVectorAD>(static_cast<std::vector<opt::ADScalar>>(cv).data(), model->nv, 1);

                    pinocchio::centerOfMass(cmodel, cdata, q_AD, false);
                    pinocchio::computeCentroidalMap(cmodel, cdata, q_AD);
                    pinocchio::forwardKinematics(cmodel, cdata, q_AD, v_AD);
                    pinocchio::updateFramePlacements(cmodel, cdata);

                    casadi::SX ground_reaction_f;
                    casadi::SX ground_reaction_tau;

                    std::vector<casadi::SX> ground_reaction_f_vec;
                    std::vector<casadi::SX> ground_reaction_tau_vec;

                    // First we create a dynamics function given only the TOTAL WRENCH

                    casadi::SX total_f(3, 1);
                    casadi::SX total_tau(3, 1);

                    auto mass = data.mass[0];
                    auto g = casadi::SX::zeros(3, 1);
                    g(2) = 9.81;
                    // derivative of moment
                    auto cdh = si->get_cdh(cx);
                    // joint velocity inputs
                    auto cvju = si->get_vju(cu);

                    // Centroidal momentum matrix
                    auto Ag = cdata.Ag;
                    casadi::SX cAg(Ag.rows(), Ag.cols());
                    pinocchio::casadi::copy(Ag, cAg);
                    auto ch = si->get_ch(cx);

                    casadi::SX cf = casadi::SX::sym("wrench_f", 3);
                    casadi::SX ctau = casadi::SX::sym("wrench_tau", 3);

                    // cu is the wrench, appended onto the joint velocities.
                    casadi::Function F("F",
                                       {cx, cf, ctau, cvju},
                                       {vertcat(cdh,
                                                (cf - mass * g) / mass,
                                                ctau / mass,
                                                cv,
                                                casadi::SX::mtimes(
                                                    casadi::SX::inv(cAg(casadi::Slice(0, 6), casadi::Slice(0, 6))),
                                                    (mass * ch - casadi::SX::mtimes(cAg(casadi::Slice(0, 6), casadi::Slice(6, int(Ag.cols()))),
                                                                                    cvju))),
                                                cvju)});

                    // Now, we create a function to gain the total wrench
                    casadi::SX total_wrench = casadi::SX::zeros(6, 1);

                    std::vector<casadi::SX> foot_forces;
                    std::vector<casadi::SX> foot_poss;
                    std::vector<casadi::SX> foot_taus;

                    for (auto ee : combination_definition)
                    {
                        if (ee.second)
                        {
                            // get the wrench "contribution" if in contact
                            auto end_effector_ptr = end_effectors[ee.first];
                            auto foot_pos = cdata.oMf[end_effector_ptr->frame_id].translation() - cdata.com[0];

                            casadi::SX cfoot_pos(3, 1);
                            pinocchio::casadi::copy(foot_pos, cfoot_pos);

                            auto cforce = si->get_f(cu, ee.first);

                            casadi::SX ctau_ee;
                            if (end_effector_ptr->is_6d)
                            {
                                ctau_ee = si->get_tau(cu, ee.first);
                            }
                            else
                            {
                                ctau_ee = casadi::SX::zeros(3, 1);
                            }

                            foot_forces.push_back(cforce);
                            foot_poss.push_back(cfoot_pos);
                            foot_taus.push_back(ctau_ee);
                        }
                    }

                    // Now we have the foot forces, foot positions, and foot torques.
                    // We can now create the total wrench, and our mode specific dynamics function

                    casadi::SX total_f_input = casadi::SX::zeros(3, 1);
                    casadi::SX total_tau_input = casadi::SX::zeros(3, 1);
                    for (std::size_t i = 0; i < foot_forces.size(); i++)
                    {
                        total_f_input += foot_forces[i];
                        total_tau_input += cross(foot_forces[i], foot_poss[i]) + foot_taus[i];
                    }

                    casadi::Function F_mode = casadi::Function("F_mode",
                                                               {cx, cu},
                                                               {F(casadi::SXVector{cx,
                                                                                   total_f_input,
                                                                                   total_tau_input,
                                                                                   cvju})
                                                                    .at(0)});
                    std::cout << "Mode dynamics created" << std::endl;
                    ModeDynamics_ = std::make_shared<casadi::Function>(F_mode);
                }

                /**
                 * @brief Gets the dynamics for this mode
                 *
                 * @param f The dynamics function
                 */
                void getModeDynamics(casadi::Function &f)
                {
                    f = *ModeDynamics_;
                }

            protected:
                std::shared_ptr<casadi::Function> ModeDynamics_;
            };

            /**
             * @brief Class for holding simple contact sequence metadata.
             *
             */
            class ContactSequence : public opt::PhaseSequence<ContactMode>
            {
            public:
                /**
                 * @brief Error codes.
                 *
                 */
                typedef PHASE_SEQUENCE_ERROR CONTACT_SEQUENCE_ERROR;

                /**
                 * @brief Construct a new Contact Sequence object.
                 *
                 * @param num_end_effectors The number of end effectors in the contact sequence.
                 */
                ContactSequence(int num_end_effectors) : PhaseSequence(), num_end_effectors_(num_end_effectors) {}

                ~ContactSequence() {}

                const int &num_end_effectors() { return num_end_effectors_; }

                /**
                 * @brief Adds a phase to the contact sequence.
                 *
                 * @param mode The contact mode of the phase.
                 * @param knot_points The number of knot points in the phase.
                 * @param dt The time step of the phase.
                 */
                int addPhase(const ContactMode &mode, int knot_points, double dt);

            private:
                /**
                 * @brief The number of end effectors in the contact sequence.
                 */
                int num_end_effectors_;
            };
        }
    }
}