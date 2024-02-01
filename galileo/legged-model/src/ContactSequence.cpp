#include "galileo/legged-model/ContactSequence.h"

namespace galileo
{
    namespace legged
    {
        namespace contact
        {
            const environment::SurfaceID &ContactMode::getSurfaceID(const std::string &ee_name) const
            {
                const auto &map_element = combination_definition.find(ee_name);

                if (map_element == combination_definition.end())
                {
                    throw std::runtime_error(std::string("'End Effector " + ee_name + " not defined in contact mode!'"));
                }

                uint index = std::distance(combination_definition.begin(), map_element);

                return contact_surfaces[index];
            }

            // Makes the combination valid. If an EE is not in contact, it makes the corresponding contact surface NO_SURFACE
            void ContactMode::MakeValid(ContactModeValidity &validity)
            {
                if (combination_definition.size() != contact_surfaces.size())
                {
                    validity = ContactModeValidity::DIFFERING_SIZES;
                    return;
                }

                ContactCombination::iterator it = combination_definition.begin();

                for (uint i = 0; i < combination_definition.size(); i++)
                {
                    bool is_in_contact = (*it).second;

                    if (!is_in_contact)
                    {
                        // If it is not in contact, make the contact surface NO_SURFACE
                        contact_surfaces[i] = environment::NO_SURFACE;
                    }

                    bool surface_is_not_defined = (contact_surfaces[i] == environment::NO_SURFACE);

                    if (is_in_contact && surface_is_not_defined)
                    {
                        validity = ContactModeValidity::SURFACE_NOT_DEFINED;
                        return;
                    }

                    if (i < combination_definition.size() - 1)
                        it++;
                }

                validity = ContactMode::ContactModeValidity::VALID;
            }

            int ContactSequence::addPhase(const ContactMode &mode, int knot_points, double dt)
            {
                std::cout << "Adding phase from Contact Sequence class" << std::endl;
                auto validity_mode = mode;
                ContactMode::ContactModeValidity validity;

                validity_mode.MakeValid(validity);

                if (validity != ContactMode::VALID)
                {
                    throw std::runtime_error(std::string("Contact is not valid!"));
                }
                std::cout << "Contact is valid!" << std::endl;
                return commonAddPhase(mode, knot_points, dt);
            }

            // void ContactMode::createModeDynamics(std::shared_ptr<opt::Model> model, RobotEndEffectors end_effectors, std::shared_ptr<opt::LeggedRobotStates> states)
            // {

            //     // opt::Data data = opt::Data(*model);

            //     // auto si = states;
            //     // casadi::SX cx = casadi::SX::sym("x", si->nx);
            //     // casadi::SX cu = casadi::SX::sym("u", si->nu);

            //     // // Get the slices of the state that represent
            //     // // Joint pos/vel
            //     // auto cq = si->get_q(cx);
            //     // auto cv = si->get_v(cx);

            //     // opt::ADModel cmodel = (*model).cast<opt::ADScalar>();
            //     // opt::ADData cdata(cmodel);
            //     // opt::ConfigVectorAD q_AD(model->nq);
            //     // opt::TangentVectorAD v_AD(model->nv);

            //     // q_AD = Eigen::Map<opt::ConfigVectorAD>(static_cast<std::vector<opt::ADScalar>>(cq).data(), model->nq, 1);
            //     // v_AD = Eigen::Map<opt::TangentVectorAD>(static_cast<std::vector<opt::ADScalar>>(cv).data(), model->nv, 1);

            //     // pinocchio::centerOfMass(cmodel, cdata, q_AD, false);
            //     // pinocchio::computeCentroidalMap(cmodel, cdata, q_AD);
            //     // pinocchio::forwardKinematics(cmodel, cdata, q_AD, v_AD);
            //     // pinocchio::updateFramePlacements(cmodel, cdata);

            //     // casadi::SX ground_reaction_f;
            //     // casadi::SX ground_reaction_tau;

            //     // std::vector<casadi::SX> ground_reaction_f_vec;
            //     // std::vector<casadi::SX> ground_reaction_tau_vec;

            //     // // First we create a dynamics function given only the TOTAL WRENCH

            //     // casadi::SX total_f(3, 1);
            //     // casadi::SX total_tau(3, 1);

            //     // auto mass = data.mass[0];
            //     // auto g = casadi::SX::zeros(3, 1);
            //     // g(2) = 9.81;
            //     // // derivative of moment
            //     // auto cdh = si->get_cdh(cx);
            //     // // joint velocity inputs
            //     // auto cvju = si->get_vju(cu);

            //     // // Centroidal momentum matrix
            //     // auto Ag = cdata.Ag;
            //     // casadi::SX cAg(Ag.rows(), Ag.cols());
            //     // pinocchio::casadi::copy(Ag, cAg);
            //     // auto ch = si->get_ch(cx);

            //     // casadi::SX cf = casadi::SX::sym("wrench_f", 3);
            //     // casadi::SX ctau = casadi::SX::sym("wrench_tau", 3);

            //     // // cu is the wrench, appended onto the joint velocities.
            //     // casadi::Function F("F",
            //     //                    {cx, cf, ctau, cvju},
            //     //                    {vertcat(cdh,
            //     //                             (cf - mass * g) / mass,
            //     //                             ctau / mass,
            //     //                             cv,
            //     //                             casadi::SX::mtimes(
            //     //                                 casadi::SX::inv(cAg(casadi::Slice(0, 6), casadi::Slice(0, 6))),
            //     //                                 (mass * ch - casadi::SX::mtimes(cAg(casadi::Slice(0, 6), casadi::Slice(6, int(Ag.cols()))),
            //     //                                                                 cvju))),
            //     //                             cvju)});

            //     // // Now, we create a function to gain the total wrench
            //     // casadi::SX total_wrench = casadi::SX::zeros(6, 1);

            //     // std::vector<casadi::SX> foot_forces;
            //     // std::vector<casadi::SX> foot_poss;
            //     // std::vector<casadi::SX> foot_taus;

            //     // for (auto ee : combination_definition)
            //     // {
            //     //     if (ee.second)
            //     //     {
            //     //         // The index of the end effector in the combination definiton
            //     //         int end_effector_idx = std::distance(combination_definition.begin(), combination_definition.find(ee.first));

            //     //         // get the wrench "contribution" if in contact
            //     //         auto end_effector_ptr = end_effectors[ee.first];
            //     //         auto foot_pos = cdata.oMf[end_effector_ptr->frame_id].translation() - cdata.com[0];

            //     //         casadi::SX cfoot_pos(3, 1);
            //     //         pinocchio::casadi::copy(foot_pos, cfoot_pos);

            //     //         auto cforce = si->get_f(cu, end_effector_idx);

            //     //         casadi::SX ctau_ee;
            //     //         if (end_effector_ptr->is_6d)
            //     //         {
            //     //             ctau_ee = si->get_tau(cu, end_effector_idx);
            //     //         }
            //     //         else
            //     //         {
            //     //             ctau_ee = casadi::SX::zeros(3, 1);
            //     //         }

            //     //         foot_forces.push_back(cforce);
            //     //         foot_poss.push_back(cfoot_pos);
            //     //         foot_taus.push_back(ctau_ee);
            //     //     }
            //     // }

            //     // // Now we have the foot forces, foot positions, and foot torques.
            //     // // We can now create the total wrench, and our mode specific dynamics function

            //     // casadi::SX total_f_input = casadi::SX::zeros(3, 1);
            //     // casadi::SX total_tau_input = casadi::SX::zeros(3, 1);
            //     // for (std::size_t i = 0; i < foot_forces.size(); i++)
            //     // {
            //     //     total_f_input += foot_forces[i];
            //     //     total_tau_input += cross(foot_forces[i], foot_poss[i]) + foot_taus[i];
            //     // }

            //     // casadi::Function F_mode = casadi::Function("F_mode",
            //     //                                            {cx, cu},
            //     //                                            {F(casadi::SXVector{cx,
            //     //                                                                total_f_input,
            //     //                                                                total_tau_input,
            //     //                                                                cvju})
            //     //                                                 .at(0)});

            //     // ModeDynamics_ = F_mode;
            // }
        }
    }
}