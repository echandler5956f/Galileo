#include "galileo/legged-model/LeggedBody.h"

namespace galileo
{
    namespace legged
    {
        LeggedBody::LeggedBody(const std::string location, const int num_ees, const std::string end_effector_names[])
        {
            model = opt::Model();
            std::vector<std::string> ee_name_vect;
            ee_name_vect.resize(num_ees);
            for (int i = 0; i < num_ees; ++i)
            {
                ee_name_vect[i] = end_effector_names[i];
            }
            pinocchio::urdf::buildModel(location, pinocchio::JointModelFreeFlyer(), model);
            data = opt::Data(model);
            cmodel = model.cast<opt::ADScalar>();
            cdata = opt::ADData(cmodel);
            setEndEffectors(ee_name_vect);
            generateContactCombination();
            si = std::make_shared<opt::LeggedRobotStates>(model.nq, model.nv, ees_);
            contact_sequence = std::make_shared<contact::ContactSequence>(num_ees);

            cx = casadi::SX::sym("x", si->nx);
            cdx = casadi::SX::sym("dx", si->ndx);
            cu = casadi::SX::sym("u", si->nu);
            cu_general = casadi::SX::sym("u_general", si->nu_general);
            cdt = casadi::SX::sym("dt");

            q_AD = Eigen::Map<opt::ConfigVectorAD>(static_cast<std::vector<opt::ADScalar>>(si->get_q(cx)).data(), model.nq, 1);
            v_AD = Eigen::Map<opt::TangentVectorAD>(static_cast<std::vector<opt::ADScalar>>(si->get_v(cx)).data(), model.nv, 1);
            createGeneralFunctions();
        }

        std::vector<pinocchio::JointIndex> getJointIndicesFromRootToFrame(const pinocchio::Model &model, pinocchio::FrameIndex frame_id)
        {
            std::vector<pinocchio::JointIndex> joint_indices;

            // Get the joint of the frame
            pinocchio::JointIndex joint_id = model.frames[frame_id].parentJoint;

            // Go up the tree until we reach the root
            while (model.parents[joint_id] > 0)
            {
                // Add the joint to the list
                joint_indices.push_back(joint_id);

                // Go to the parent joint
                joint_id = model.parents[joint_id];
            }

            // Reverse the list to get the joints from the root to the frame
            std::reverse(joint_indices.begin(), joint_indices.end());

            return joint_indices;
        }

        void LeggedBody::setEndEffectors(const std::vector<std::string> &ee_names)
        {
            std::vector<pinocchio::FrameIndex> ee_ids;
            for (std::size_t i = 0; i < ee_names.size(); ++i)
            {
                auto ee_name = ee_names[i];
                assert(model.existFrame(ee_name, pinocchio::BODY));
                // Throw an error otherwise.
                std::shared_ptr<contact::EndEffector> ee_obj_ptr = std::make_shared<contact::EndEffector>();
                ee_obj_ptr->frame_name = ee_name;
                pinocchio::FrameIndex frame_id = model.getFrameId(ee_name);
                ee_obj_ptr->frame_id = frame_id;

                // Get the joint indices from the root to the end effector
                std::vector<pinocchio::JointIndex> joint_indices = getJointIndicesFromRootToFrame(model, frame_id);

                // Compute the total number of degrees of freedom
                int dof = 0;
                for (const auto &joint_index : joint_indices)
                {
                    dof += model.joints[joint_index].nq();
                }
                ee_obj_ptr->is_6d = dof == 6;
                ee_obj_ptr->local_ee_idx = i;

                ees_.insert({frame_id, ee_obj_ptr});
                ee_ids.push_back(frame_id);
            }
            num_end_effectors_ = ee_names.size();
            ee_ids_ = ee_ids;
        }

        void LeggedBody::generateContactCombination()
        {
            std::vector<contact::ContactCombination> contact_combinations;
            // Generate the "basic" (no contact) contact combination.
            contact::ContactCombination basic_cc;
            for (auto &ee_id : ee_ids_)
                basic_cc.insert({ee_id, false});

            // The power set can be gotten by all the binary values between 0 and 2^n-1.
            // That is, 0000, 0001, 0010, 0011, ... for 4 EEs.
            uint num_combinations = pow(2, num_end_effectors_);

            contact_combinations.resize(num_combinations);
            for (uint binary_value_combination = 0; binary_value_combination < num_combinations; ++binary_value_combination)
            {
                // Copy the no contact cc
                contact::ContactCombination new_contact_combination = basic_cc;
                // And set the value for each ee in contact to true
                std::bitset<32> bit_mask = 1;
                std::bitset<32> binary_value_combination_bits = binary_value_combination;

                for (int i = 0; i < num_end_effectors_; i++)
                {
                    bool ee_i_is_in_contact = (bit_mask & binary_value_combination_bits).any();
                    // bit shift
                    bit_mask <<= 1;
                    new_contact_combination[ee_ids_[i]] = ee_i_is_in_contact;
                }
                contact_combinations[binary_value_combination] = new_contact_combination;
            }
            contact_combinations_ = contact_combinations;
        }

        void LeggedBody::createGeneralFunctions()
        {
            createGeneralDynamics();
            createFint();
            createFdiff();
        }

        void LeggedBody::createGeneralDynamics()
        {
            pinocchio::centerOfMass(cmodel, cdata, q_AD, false);
            pinocchio::computeCentroidalMap(cmodel, cdata, q_AD);
            pinocchio::forwardKinematics(cmodel, cdata, q_AD, v_AD);
            pinocchio::updateFramePlacements(cmodel, cdata);

            casadi::SX mass = cdata.mass[0];
            casadi::SX g = casadi::SX::zeros(3, 1);
            g(2) = 9.81;

            auto Ag = cdata.Ag;
            casadi::SX cAg(Ag.rows(), Ag.cols());
            pinocchio::casadi::copy(Ag, cAg);

            general_dynamics = casadi::Function("F",
                                                {cx, cu_general},
                                                {vertcat(si->get_cdh(cx),
                                                         (si->get_general_forces(cu_general) - mass * g) / mass,
                                                         si->get_general_torques(cu_general) / mass,
                                                         si->get_v(cx),
                                                         casadi::SX::mtimes(
                                                             casadi::SX::inv(cAg(casadi::Slice(0, cAg.size1()), casadi::Slice(0, 6))),
                                                             (mass * si->get_ch(cx) - casadi::SX::mtimes(cAg(casadi::Slice(0, cAg.size1()), casadi::Slice(6, cAg.size2())),
                                                                                                         si->get_general_joint_velocities(cu_general)))),
                                                         si->get_general_joint_velocities(cu_general))});
        }

        void LeggedBody::createFint()
        {
            casadi::SX qb = si->get_q(cx)(casadi::Slice(0, si->nqb));
            casadi::SX vb = si->get_q_d(cdx)(casadi::Slice(0, si->nvb));
            casadi::SX lie_group_int_result = lie_group_int(qb, vb, cdt);
            casadi::SX q_joints_int_result = si->get_q(cx)(casadi::Slice(si->nqb, si->nq)) + si->get_q_d(cdx)(casadi::Slice(si->nvb, si->nv)) * cdt;

            fint = casadi::Function("Fint",
                                    {cx, cdx, cdt},
                                    {vertcat(si->get_ch(cx) + si->get_ch_d(cdx) * cdt,
                                             si->get_cdh(cx) + si->get_cdh_d(cdx) * cdt,
                                             lie_group_int_result,
                                             q_joints_int_result,
                                             si->get_v(cx) + si->get_v_d(cdx) * cdt)});
        }

        casadi::SX LeggedBody::lie_group_int(casadi::SX qb, casadi::SX vb, casadi::SX dt)
        {
            casadi::SX pos = qb(casadi::Slice(0, 3));
            casadi::SX quat = qb(casadi::Slice(3, 7));

            casadi::SX vel = vb(casadi::Slice(0, 3)) * dt;
            casadi::SX omega = vb(casadi::Slice(3, 6)) * dt;
            casadi::SX omega_quat = casadi::SX::vertcat({omega / 2, 0});
            casadi::SX exp_omega_quat = quat_exp(omega_quat);
            casadi::SX quat_new = quat_mult(quat, exp_omega_quat);
            quat_new = quat_new / sqrt(quat_new(0) * quat_new(0) + quat_new(1) * quat_new(1) + quat_new(2) * quat_new(2) + quat_new(3) * quat_new(3) + casadi::eps * casadi::eps);
            casadi::SX pos_new = pos + apply_quat(quat, casadi::SX::mtimes(rodrigues(omega), vel));
            return vertcat(pos_new, quat_new);
        }

        casadi::SX LeggedBody::quat_exp(casadi::SX quat)
        {
            casadi::SX v = quat(casadi::Slice(0, 3));
            casadi::SX v_norm = sqrt(v(0) * v(0) + v(1) * v(1) + v(2) * v(2) + casadi::eps * casadi::eps);
            casadi::SX a = quat(3);

            casadi::SX exp_a = exp(a);
            casadi::SX cos_norm_v = cos(v_norm);
            casadi::SX sin_norm_v = sin(v_norm);

            casadi::SX quat_new = casadi::SX::vertcat({exp_a * v * sin_norm_v / v_norm, exp_a * cos_norm_v});

            return quat_new / sqrt(quat_new(0) * quat_new(0) + quat_new(1) * quat_new(1) + quat_new(2) * quat_new(2) + quat_new(3) * quat_new(3) + casadi::eps * casadi::eps);
        }

        void LeggedBody::createFdiff()
        {
            casadi::SX qb = si->get_q(cx)(casadi::Slice(0, si->nqb));
            casadi::SX cx2 = casadi::SX::sym("x2", si->nx);
            casadi::SX qb2 = si->get_q(cx2)(casadi::Slice(0, si->nqb));
            casadi::SX lie_group_diff_result = lie_group_diff(qb, qb2, cdt);
            casadi::SX v_joints_diff_result = (si->get_vj(cx2) - si->get_vj(cx)) / cdt;

            fdiff = casadi::Function("Fdiff",
                                     {cx, cx2, cdt},
                                     {vertcat((si->get_ch(cx2) - si->get_ch(cx)) / cdt,
                                              (si->get_cdh(cx2) - si->get_cdh(cx)) / cdt,
                                              lie_group_diff_result,
                                              v_joints_diff_result,
                                              (si->get_v(cx2) - si->get_v(cx)) / cdt)});
        }

        casadi::SX LeggedBody::lie_group_diff(casadi::SX qb1, casadi::SX qb2, casadi::SX dt)
        {
            casadi::SX pos1 = qb1(casadi::Slice(0, 3));
            casadi::SX quat1 = qb1(casadi::Slice(3, 7));
            casadi::SX pos2 = qb2(casadi::Slice(0, 3));
            casadi::SX quat2 = qb2(casadi::Slice(3, 7));

            casadi::SX quat_diff = quat_mult(quat_inv(quat1), quat2);
            casadi::SX omega_quat = quat_log(quat_diff);
            casadi::SX omega = omega_quat(casadi::Slice(0, 3));

            casadi::SX pos_diff = pos2 - pos1;
            casadi::SX vel = casadi::SX::mtimes(inv(rodrigues(omega)), apply_quat(quat_inv(quat1), pos_diff));

            return casadi::SX::vertcat({vel / dt, omega / dt});
        }

        casadi::SX LeggedBody::quat_log(casadi::SX quat)
        {
            casadi::SX squared_n = quat(0) * quat(0) + quat(1) * quat(1) + quat(2) * quat(2);
            casadi::SX n = sqrt(squared_n);
            casadi::SX w = quat(3);

            return casadi::SX::if_else(n < casadi::eps,
                                                ((2./w) - (2./3.) * (squared_n / (w * w * w))) * quat(casadi::Slice(0, 3)),
                                                4. * atan(n / (w + sqrt(w * w + squared_n + casadi::eps * casadi::eps))) * quat(casadi::Slice(0, 3)) / n);
        }

        casadi::SX LeggedBody::apply_quat(casadi::SX quat, casadi::SX vec3)
        {
            casadi::SX imag = quat(casadi::Slice(0, 3));
            casadi::SX real = quat(3);
            return vec3 + 2 * cross(imag, cross(imag, vec3) + real * vec3);
        }

        casadi::SX LeggedBody::quat_mult(casadi::SX quat1, casadi::SX quat2)
        {
            casadi::SX real_1 = quat1(3);
            casadi::SX imag_1 = quat1(casadi::Slice(0, 3));

            casadi::SX real_2 = quat2(3);
            casadi::SX imag_2 = quat2(casadi::Slice(0, 3));

            casadi::SX real_res = real_1 * real_2 - dot(imag_1, imag_2);
            casadi::SX imag_res = real_1 * imag_2 + real_2 * imag_1 + cross(imag_1, imag_2);
            return casadi::SX::vertcat({imag_res, real_res});
        }

        casadi::SX LeggedBody::quat_inv(casadi::SX quat)
        {
            casadi::SX quat_norm = sqrt(quat(0) * quat(0) + quat(1) * quat(1) + quat(2) * quat(2) + quat(3) * quat(3) + casadi::eps * casadi::eps);
            return casadi::SX::vertcat({-quat(casadi::Slice(0, 3)) / quat_norm, quat(3) / quat_norm});
        }

        casadi::SX LeggedBody::rodrigues(casadi::SX omega)
        {
            casadi::SX theta = sqrt(omega(0) * omega(0) + omega(1) * omega(1) + omega(2) * omega(2) + casadi::eps * casadi::eps);
            return casadi::SX::eye(3) + ((1 - cos(theta)) / (theta * theta)) * skew(omega) + ((theta - sin(theta)) / (theta * theta * theta)) * mpower(skew(omega), 2);
        }

        void LeggedBody::fillModeDynamics(bool print_ees_info)
        {
            casadi::SXVector foot_forces;
            casadi::SXVector foot_poss;
            casadi::SXVector foot_taus;
            for (std::size_t i = 0; i < contact_sequence->phase_sequence_.size(); ++i)
            {
                contact::ContactMode mode = contact_sequence->phase_sequence_[i].mode;

                foot_forces.clear();
                foot_poss.clear();
                foot_taus.clear();

                for (auto ee : mode.combination_definition)
                {
                    std::string print_info;
                    auto end_effector_ptr = ees_[ee.first];
                    if (ee.second)
                    {
                        print_info = "At phase " + std::to_string(i) + " " + end_effector_ptr->frame_name + " is in contact.";
                        auto foot_pos = cdata.oMf[end_effector_ptr->frame_id].translation() - cdata.com[0];

                        casadi::SX cfoot_pos(3, 1);
                        pinocchio::casadi::copy(foot_pos, cfoot_pos);

                        casadi::SX cforce = si->get_f(cu, ee.first);

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
                    else
                        print_info = "At phase " + std::to_string(i) + " " + end_effector_ptr->frame_name + " is NOT in contact.";

                    if (print_ees_info)
                        std::cout << print_info << std::endl;
                }

                casadi::SX total_f_input = casadi::SX::zeros(3, 1);
                casadi::SX total_tau_input = casadi::SX::zeros(3, 1);
                for (std::size_t j = 0; j < foot_forces.size(); j++)
                {
                    total_f_input += foot_forces[j];
                    total_tau_input += cross(foot_poss[j], foot_forces[j]) + foot_taus[j];
                }

                casadi::SX u_general = vertcat(casadi::SXVector{total_f_input, total_tau_input, si->get_vju(cu)});

                contact_sequence->phase_sequence_[i].phase_dynamics = casadi::Function("F_mode",
                                                                                       {cx, cu},
                                                                                       {general_dynamics(casadi::SXVector{cx, u_general})
                                                                                            .at(0)});
            }
        }

        contact::ContactCombination LeggedBody::getContactCombination(int contact_mask) { return contact_combinations_[contact_mask]; }

        contact::RobotEndEffectors LeggedBody::getEndEffectors() { return ees_; }
    }
}
