#include "galileo/legged-model/LeggedBody.h"

namespace galileo
{
    namespace legged
    {
        LeggedBody::LeggedBody(const std::string location, const std::vector<std::string> end_effector_names, casadi::Dict general_function_casadi_options)
        {
            model = Model();

            pinocchio::urdf::buildModel(location, pinocchio::JointModelFreeFlyer(), model);
            data = Data(model);
            cmodel = model.cast<ADScalar>();
            cdata = ADData(cmodel);
            setEndEffectors(end_effector_names);
            generateContactCombination();
            si = std::make_shared<legged::LeggedRobotStates>(model.nq, model.nv, ees_);
            contact_sequence = std::make_shared<contact::ContactSequence>(end_effector_names.size());

            cx = casadi::SX::sym("x", si->nx);
            cdx = casadi::SX::sym("dx", si->ndx);
            cu = casadi::SX::sym("u", si->nu);
            cu_general = casadi::SX::sym("u_general", si->nu_general);
            cdt = casadi::SX::sym("dt");

            q_AD = Eigen::Map<ConfigVectorAD>(static_cast<std::vector<ADScalar>>(si->get_q(cx)).data(), model.nq, 1);
            pinocchio::forwardKinematics(cmodel, cdata, q_AD);
            pinocchio::framesForwardKinematics(cmodel, cdata, q_AD);
            pinocchio::updateFramePlacements(cmodel, cdata);

            createGeneralFunctions(general_function_casadi_options);
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

        Eigen::MatrixXd LeggedBody::initializeInputCostWeight(Eigen::MatrixXd R_taskspace, ConfigVector q0)
        {
            pinocchio::computeJointJacobians(model, data, q0);
            pinocchio::updateFramePlacements(model, data);

            size_t totalContactDim = 3 * num_end_effectors_;

            Eigen::MatrixXd baseToFeetJacobians(totalContactDim, si->nvju);
            int i = 0;
            for (auto ee : ees_)
            {
                auto frame_id = ee.first;
                Eigen::MatrixXd jacobianWorldToContactPointInWorldFrame = Eigen::MatrixXd::Zero(6, si->nv);
                pinocchio::getFrameJacobian(model, data, frame_id, pinocchio::LOCAL_WORLD_ALIGNED, jacobianWorldToContactPointInWorldFrame);
                baseToFeetJacobians.block(3 * i, 0, 3, si->nvju) =
                    jacobianWorldToContactPointInWorldFrame.block(0, 6, 3, si->nvju);
                ++i;
            }

            Eigen::MatrixXd R = Eigen::MatrixXd::Zero(si->nu, si->nu);
            // Contact Forces
            R.topLeftCorner(totalContactDim, totalContactDim) = R_taskspace.topLeftCorner(totalContactDim, totalContactDim);
            // Joint velocities
            R.bottomRightCorner(si->nvju, si->nvju) = baseToFeetJacobians.transpose() * R_taskspace.bottomRightCorner(totalContactDim, totalContactDim) * baseToFeetJacobians;
            return R;
        }

        void LeggedBody::createGeneralFunctions(casadi::Dict casadi_opts)
        {
            createGeneralDynamics(casadi_opts);
            createFint(casadi_opts);
            createFdiff(casadi_opts);
            createErrorFunction(casadi_opts);
        }

        template <typename SCALAR_T>
        Eigen::Matrix<SCALAR_T, 6, 6> computeFloatingBaseCentroidalMomentumMatrixInverse(const Eigen::Matrix<SCALAR_T, 6, 6> &Ab)
        {
            const SCALAR_T mass = Ab(0, 0);
            Eigen::Matrix<SCALAR_T, 3, 3> Ab_22_inv = Ab.template block<3, 3>(3, 3).inverse();
            Eigen::Matrix<SCALAR_T, 6, 6> Ab_inv = Eigen::Matrix<SCALAR_T, 6, 6>::Zero();
            Ab_inv << 1.0 / mass * Eigen::Matrix<SCALAR_T, 3, 3>::Identity(), -1.0 / mass * Ab.template block<3, 3>(0, 3) * Ab_22_inv,
                Eigen::Matrix<SCALAR_T, 3, 3>::Zero(), Ab_22_inv;
            return Ab_inv;
        }

        void LeggedBody::createGeneralDynamics(casadi::Dict casadi_opts)
        {
            pinocchio::centerOfMass(cmodel, cdata, q_AD, false);
            pinocchio::computeCentroidalMap(cmodel, cdata, q_AD);

            casadi::SX mass = cdata.mass[0];
            casadi::SX g = casadi::SX::zeros(3, 1);
            g(2) = 9.81;

            TangentVectorAD vju_general_AD = Eigen::Map<TangentVectorAD>(static_cast<std::vector<ADScalar>>(si->get_general_joint_velocities(cu_general)).data(), si->nvju, 1);
            TangentVectorAD vju_AD = Eigen::Map<TangentVectorAD>(static_cast<std::vector<ADScalar>>(si->get_vju(cu)).data(), si->nvju, 1);

            auto Ag = cdata.Ag;
            const Eigen::Matrix<ADScalar, 6, 6> Ab = Ag.template leftCols<6>();
            const auto Ab_inv = computeFloatingBaseCentroidalMomentumMatrixInverse(Ab);
            TangentVectorAD h_AD = Eigen::Map<TangentVectorAD>(static_cast<std::vector<ADScalar>>(si->get_ch(cx)).data(), si->nh, 1);

            TangentVectorAD vb_AD1 = Ab_inv * (mass * h_AD - Ag.rightCols(si->nvju) * vju_general_AD);
            TangentVectorAD tmp_v_AD1(si->nv, 1);
            tmp_v_AD1 << vb_AD1, vju_general_AD;
            casadi::SX cv(si->nv, 1);
            pinocchio::casadi::copy(tmp_v_AD1, cv);

            general_dynamics = casadi::Function("F",
                                                {cx, cu_general},
                                                {vertcat((si->get_general_forces(cu_general) - mass * g) / mass,
                                                         si->get_general_torques(cu_general) / mass,
                                                         cv)});

            TangentVectorAD vb_AD2 = Ab_inv * (mass * h_AD - Ag.rightCols(si->nvju) * vju_AD);
            TangentVectorAD tmp_v_AD2(si->nv, 1);
            tmp_v_AD2 << vb_AD2, vju_AD;

            pinocchio::forwardKinematics(cmodel, cdata, q_AD, tmp_v_AD2);
            pinocchio::updateFramePlacements(cmodel, cdata);
        }

        void LeggedBody::createFint(casadi::Dict casadi_opts)
        {
            casadi::SX qb = si->get_q(cx)(casadi::Slice(0, si->nqb));
            casadi::SX vb = si->get_q_d(cdx)(casadi::Slice(0, si->nvb));
            casadi::SX lie_group_int_result = math::lie_group_int(qb, vb, cdt);
            casadi::SX q_joints_int_result = si->get_q(cx)(casadi::Slice(si->nqb, si->nq)) + si->get_q_d(cdx)(casadi::Slice(si->nvb, si->nv)) * cdt;

            fint = casadi::Function("Fint",
                                    {cx, cdx, cdt},
                                    {vertcat(si->get_ch(cx) + si->get_ch_d(cdx) * cdt,
                                             lie_group_int_result,
                                             q_joints_int_result)},
                                    casadi_opts);
        }

        void LeggedBody::createFdiff(casadi::Dict casadi_opts)
        {
            casadi::SX qb = si->get_q(cx)(casadi::Slice(0, si->nqb));
            casadi::SX cx2 = casadi::SX::sym("x2", si->nx);
            casadi::SX qb2 = si->get_q(cx2)(casadi::Slice(0, si->nqb));
            casadi::SX lie_group_diff_result = math::lie_group_diff(qb, qb2, cdt);

            fdiff = casadi::Function("Fdiff",
                                     {cx, cx2, cdt},
                                     {vertcat((si->get_ch(cx2) - si->get_ch(cx)) / cdt,
                                              lie_group_diff_result,
                                              (si->get_qj(cx2) - si->get_qj(cx)) / cdt)},
                                     casadi_opts);
        }

        void LeggedBody::createErrorFunction(casadi::Dict casadi_opts)
        {
            casadi::SX qb = si->get_q(cx)(casadi::Slice(0, si->nqb));
            casadi::SX cx2 = casadi::SX::sym("x2", si->nx);
            casadi::SX qb2 = si->get_q(cx2)(casadi::Slice(0, si->nqb));

            casadi::SX quat1 = qb(casadi::Slice(3, si->nqb));
            casadi::SX quat2 = qb2(casadi::Slice(3, si->nqb));

            casadi::SX quat_distance_result = math::quat_distance(quat1, quat2);

            f_state_error = casadi::Function("F_state_error",
                                             {cx, cx2},
                                             {vertcat(si->get_ch(cx2) - si->get_ch(cx),
                                                      qb(casadi::Slice(0, 3)) - qb2(casadi::Slice(0, 3)),
                                                      quat_distance_result,
                                                      si->get_qj(cx2) - si->get_qj(cx))},
                                             casadi_opts);
        }

        void LeggedBody::fillModeDynamics(casadi::Dict casadi_opts)
        {
            casadi::SXVector foot_forces;
            casadi::SXVector foot_poss;
            casadi::SXVector foot_taus;
            for (std::size_t i = 0; i < contact_sequence->getPhases().size(); ++i)
            {
                contact::ContactMode mode = contact_sequence->getPhases()[i].mode;

                foot_forces.clear();
                foot_poss.clear();
                foot_taus.clear();

                for (auto ee : mode.combination_definition)
                {
                    auto end_effector_ptr = ees_[ee.first];
                    if (ee.second)
                    {
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
                }

                casadi::SX total_f_input = casadi::SX::zeros(3, 1);
                casadi::SX total_tau_input = casadi::SX::zeros(3, 1);
                for (std::size_t j = 0; j < foot_forces.size(); j++)
                {
                    total_f_input += foot_forces[j];
                    total_tau_input += cross(foot_poss[j], foot_forces[j]) + foot_taus[j];
                }

                casadi::SX u_general = vertcat(casadi::SXVector{total_f_input, total_tau_input, si->get_vju(cu)});

                contact_sequence->FillPhaseDynamics(i, casadi::Function("F_mode",
                                                                        {cx, cu},
                                                                        {general_dynamics(casadi::SXVector{cx, u_general})
                                                                             .at(0)}));
            }
        }

        casadi::SX LeggedBody::weightCompensatingInputsForPhase(size_t phase_index)
        {
            casadi::SX weight_compensating_inputs = casadi::SX::zeros(si->nu, 1);
            contact::ContactMode mode = contact_sequence->getPhases()[phase_index].mode;
            for (auto ee : ees_)
            {
                if (mode[(*ee.second)])
                {
                    weight_compensating_inputs(casadi::Slice(std::get<0>(si->frame_id_to_index_range[ee.second->frame_id]) + 2)) = 9.81 * cdata.mass[0] / contact_sequence->numEndEffectorsInContactAtPhase(phase_index);
                }
            }
            return weight_compensating_inputs;
        }

        contact::ContactCombination LeggedBody::getContactCombination(int contact_mask) { return contact_combinations_[contact_mask]; }

        contact::RobotEndEffectors LeggedBody::getEndEffectors() { return ees_; }
    }
}
