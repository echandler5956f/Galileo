#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include <pinocchio/multibody/data.hpp>

#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/contact-dynamics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <pinocchio/algorithm/rnea.hpp>

#include <iostream>

using Model_t = pinocchio::ModelTpl<double, 0>;
using Data_t = pinocchio::DataTpl<double, 0>;

std::size_t get_model_size(const Model_t & model)
{
    std::size_t total_size = 0;

    // Basic int members (stack allocated)
    total_size += sizeof(model.nq);
    total_size += sizeof(model.nv);
    total_size += sizeof(model.nvExtended);
    total_size += sizeof(model.njoints);
    total_size += sizeof(model.nbodies);
    total_size += sizeof(model.nframes);

    // Vector containers - need to account for both container overhead and element storage

    // InertiaVector inertias
    total_size += sizeof(model.inertias);
    total_size += model.inertias.capacity() * sizeof(decltype(model.inertias)::value_type);

    // SE3Vector jointPlacements
    total_size += sizeof(model.jointPlacements);
    total_size += model.jointPlacements.capacity() * sizeof(decltype(model.jointPlacements)::value_type);

    // JointModelVector joints
    total_size += sizeof(model.joints);
    total_size += model.joints.capacity() * sizeof(decltype(model.joints)::value_type);

    // std::vector<int> containers
    total_size += sizeof(model.idx_qs) + model.idx_qs.capacity() * sizeof(int);
    total_size += sizeof(model.nqs) + model.nqs.capacity() * sizeof(int);
    total_size += sizeof(model.idx_vs) + model.idx_vs.capacity() * sizeof(int);
    total_size += sizeof(model.nvs) + model.nvs.capacity() * sizeof(int);
    total_size += sizeof(model.idx_vExtendeds) + model.idx_vExtendeds.capacity() * sizeof(int);
    total_size += sizeof(model.nvExtendeds) + model.nvExtendeds.capacity() * sizeof(int);

    // std::vector<JointIndex> containers
    total_size += sizeof(model.parents) + model.parents.capacity() * sizeof(pinocchio::JointIndex);
    total_size += sizeof(model.mimicking_joints) + model.mimicking_joints.capacity() * sizeof(pinocchio::JointIndex);
    total_size += sizeof(model.mimicked_joints) + model.mimicked_joints.capacity() * sizeof(pinocchio::JointIndex);

    // std::vector<IndexVector> containers - need to account for nested vectors
    total_size += sizeof(model.children);
    for (const auto& child_vec : model.children) {
        total_size += sizeof(child_vec) + child_vec.capacity() * sizeof(pinocchio::Index);
    }

    total_size += sizeof(model.supports);
    for (const auto& support_vec : model.supports) {
        total_size += sizeof(support_vec) + support_vec.capacity() * sizeof(pinocchio::Index);
    }

    total_size += sizeof(model.mimic_joint_supports);
    for (const auto& mimic_support_vec : model.mimic_joint_supports) {
        total_size += sizeof(mimic_support_vec) + mimic_support_vec.capacity() * sizeof(pinocchio::Index);
    }

    total_size += sizeof(model.subtrees);
    for (const auto& subtree_vec : model.subtrees) {
        total_size += sizeof(subtree_vec) + subtree_vec.capacity() * sizeof(pinocchio::Index);
    }

    // std::vector<std::string> names
    total_size += sizeof(model.names);
    for (const auto& name : model.names) {
        total_size += sizeof(name) + name.capacity() * sizeof(char);
    }

    // ConfigVectorMap referenceConfigurations - std::map<std::string, ConfigVectorType>
    total_size += sizeof(model.referenceConfigurations);
    for (const auto& config_pair : model.referenceConfigurations) {
        // Map node overhead (approximate)
        total_size += sizeof(std::pair<const std::string, typename Model_t::ConfigVectorType>) + 3 * sizeof(void*);
        // String key
        total_size += config_pair.first.capacity() * sizeof(char);
        // Eigen vector value
        total_size += config_pair.second.size() * sizeof(double);
    }

    // Eigen VectorXs types (dynamic vectors)
    total_size += sizeof(model.armature) + model.armature.size() * sizeof(double);
    total_size += sizeof(model.rotorInertia) + model.rotorInertia.size() * sizeof(double);
    total_size += sizeof(model.rotorGearRatio) + model.rotorGearRatio.size() * sizeof(double);
    total_size += sizeof(model.friction) + model.friction.size() * sizeof(double);
    total_size += sizeof(model.damping) + model.damping.size() * sizeof(double);
    total_size += sizeof(model.effortLimit) + model.effortLimit.size() * sizeof(double);
    total_size += sizeof(model.velocityLimit) + model.velocityLimit.size() * sizeof(double);
    total_size += sizeof(model.lowerPositionLimit) + model.lowerPositionLimit.size() * sizeof(double);
    total_size += sizeof(model.upperPositionLimit) + model.upperPositionLimit.size() * sizeof(double);

    // FrameVector frames
    total_size += sizeof(model.frames);
    total_size += model.frames.capacity() * sizeof(decltype(model.frames)::value_type);

    // Motion gravity (assuming it's a fixed-size type)
    total_size += sizeof(model.gravity);

    // std::string name
    total_size += sizeof(model.name) + model.name.capacity() * sizeof(char);

    // Note: gravity981 is static const, so it doesn't contribute to instance size

    return total_size;
}

std::size_t get_data_size(const Data_t & data)
{
    std::size_t total_size = 0;

    // JointDataVector - aligned vector of joint data
    total_size += sizeof(data.joints);
    total_size += data.joints.capacity() * sizeof(decltype(data.joints)::value_type);

    // Motion vectors (PINOCCHIO_ALIGNED_STD_VECTOR(Motion))
    total_size += sizeof(data.a) + data.a.capacity() * sizeof(decltype(data.a)::value_type);
    total_size += sizeof(data.oa) + data.oa.capacity() * sizeof(decltype(data.oa)::value_type);
    total_size += sizeof(data.oa_drift) + data.oa_drift.capacity() * sizeof(decltype(data.oa_drift)::value_type);
    total_size += sizeof(data.oa_augmented) + data.oa_augmented.capacity() * sizeof(decltype(data.oa_augmented)::value_type);
    total_size += sizeof(data.a_gf) + data.a_gf.capacity() * sizeof(decltype(data.a_gf)::value_type);
    total_size += sizeof(data.oa_gf) + data.oa_gf.capacity() * sizeof(decltype(data.oa_gf)::value_type);
    total_size += sizeof(data.v) + data.v.capacity() * sizeof(decltype(data.v)::value_type);
    total_size += sizeof(data.ov) + data.ov.capacity() * sizeof(decltype(data.ov)::value_type);
    total_size += sizeof(data.a_bias) + data.a_bias.capacity() * sizeof(decltype(data.a_bias)::value_type);

    // Force vectors (PINOCCHIO_ALIGNED_STD_VECTOR(Force))
    total_size += sizeof(data.f) + data.f.capacity() * sizeof(decltype(data.f)::value_type);
    total_size += sizeof(data.of) + data.of.capacity() * sizeof(decltype(data.of)::value_type);
    total_size += sizeof(data.of_augmented) + data.of_augmented.capacity() * sizeof(decltype(data.of_augmented)::value_type);
    total_size += sizeof(data.h) + data.h.capacity() * sizeof(decltype(data.h)::value_type);
    total_size += sizeof(data.oh) + data.oh.capacity() * sizeof(decltype(data.oh)::value_type);

    // SE3 vectors (PINOCCHIO_ALIGNED_STD_VECTOR(SE3))
    total_size += sizeof(data.oMi) + data.oMi.capacity() * sizeof(decltype(data.oMi)::value_type);
    total_size += sizeof(data.liMi) + data.liMi.capacity() * sizeof(decltype(data.liMi)::value_type);
    total_size += sizeof(data.oMf) + data.oMf.capacity() * sizeof(decltype(data.oMf)::value_type);
    total_size += sizeof(data.iMf) + data.iMf.capacity() * sizeof(decltype(data.iMf)::value_type);

    // Eigen VectorXs and TangentVectorType (dynamic vectors)
    total_size += sizeof(data.tau) + data.tau.size() * sizeof(double);
    total_size += sizeof(data.nle) + data.nle.size() * sizeof(double);
    total_size += sizeof(data.g) + data.g.size() * sizeof(double);
    total_size += sizeof(data.ddq) + data.ddq.size() * sizeof(double);
    total_size += sizeof(data.u) + data.u.size() * sizeof(double);
    total_size += sizeof(data.D) + data.D.size() * sizeof(double);
    total_size += sizeof(data.Dinv) + data.Dinv.size() * sizeof(double);
    total_size += sizeof(data.tmp) + data.tmp.size() * sizeof(double);
    total_size += sizeof(data.lambda_c) + data.lambda_c.size() * sizeof(double);
    total_size += sizeof(data.lambda_c_prox) + data.lambda_c_prox.size() * sizeof(double);
    total_size += sizeof(data.diff_lambda_c) + data.diff_lambda_c.size() * sizeof(double);
    total_size += sizeof(data.torque_residual) + data.torque_residual.size() * sizeof(double);
    total_size += sizeof(data.dq_after) + data.dq_after.size() * sizeof(double);
    total_size += sizeof(data.impulse_c) + data.impulse_c.size() * sizeof(double);
    total_size += sizeof(data.primal_dual_contact_solution) + data.primal_dual_contact_solution.size() * sizeof(double);
    total_size += sizeof(data.primal_rhs_contact) + data.primal_rhs_contact.size() * sizeof(double);

    // Inertia vectors (PINOCCHIO_ALIGNED_STD_VECTOR(Inertia))
    total_size += sizeof(data.Ycrb) + data.Ycrb.capacity() * sizeof(decltype(data.Ycrb)::value_type);
    total_size += sizeof(data.oinertias) + data.oinertias.capacity() * sizeof(decltype(data.oinertias)::value_type);
    total_size += sizeof(data.oYcrb) + data.oYcrb.capacity() * sizeof(decltype(data.oYcrb)::value_type);

    // Matrix6 vectors (PINOCCHIO_ALIGNED_STD_VECTOR(Matrix6))
    total_size += sizeof(data.dYcrb) + data.dYcrb.capacity() * sizeof(decltype(data.dYcrb)::value_type);
    total_size += sizeof(data.vxI) + data.vxI.capacity() * sizeof(decltype(data.vxI)::value_type);
    total_size += sizeof(data.Ivx) + data.Ivx.capacity() * sizeof(decltype(data.Ivx)::value_type);
    total_size += sizeof(data.B) + data.B.capacity() * sizeof(decltype(data.B)::value_type);
    total_size += sizeof(data.doYcrb) + data.doYcrb.capacity() * sizeof(decltype(data.doYcrb)::value_type);
    total_size += sizeof(data.Yaba) + data.Yaba.capacity() * sizeof(decltype(data.Yaba)::value_type);
    total_size += sizeof(data.oYaba) + data.oYaba.capacity() * sizeof(decltype(data.oYaba)::value_type);
    total_size += sizeof(data.oYaba_contact) + data.oYaba_contact.capacity() * sizeof(decltype(data.oYaba_contact)::value_type);
    total_size += sizeof(data.oL) + data.oL.capacity() * sizeof(decltype(data.oL)::value_type);
    total_size += sizeof(data.oK) + data.oK.capacity() * sizeof(decltype(data.oK)::value_type);
    total_size += sizeof(data.extended_motion_propagator) + data.extended_motion_propagator.capacity() * sizeof(decltype(data.extended_motion_propagator)::value_type);
    total_size += sizeof(data.extended_motion_propagator2) + data.extended_motion_propagator2.capacity() * sizeof(decltype(data.extended_motion_propagator2)::value_type);
    total_size += sizeof(data.spatial_inv_inertia) + data.spatial_inv_inertia.capacity() * sizeof(decltype(data.spatial_inv_inertia)::value_type);

    // Eigen MatrixXs (dynamic matrices)
    total_size += sizeof(data.M) + data.M.size() * sizeof(double);
    total_size += sizeof(data.Minv) + data.Minv.size() * sizeof(double);
    total_size += sizeof(data.C) + data.C.size() * sizeof(double);
    total_size += sizeof(data.U) + data.U.size() * sizeof(double);
    total_size += sizeof(data.JMinvJt) + data.JMinvJt.size() * sizeof(double);
    total_size += sizeof(data.sDUiJt) + data.sDUiJt.size() * sizeof(double);
    total_size += sizeof(data.dvc_dq) + data.dvc_dq.size() * sizeof(double);
    total_size += sizeof(data.dac_dq) + data.dac_dq.size() * sizeof(double);
    total_size += sizeof(data.dac_dv) + data.dac_dv.size() * sizeof(double);
    total_size += sizeof(data.dac_da) + data.dac_da.size() * sizeof(double);
    total_size += sizeof(data.osim) + data.osim.size() * sizeof(double);
    total_size += sizeof(data.dlambda_dq) + data.dlambda_dq.size() * sizeof(double);
    total_size += sizeof(data.dlambda_dv) + data.dlambda_dv.size() * sizeof(double);
    total_size += sizeof(data.dlambda_dtau) + data.dlambda_dtau.size() * sizeof(double);
    total_size += sizeof(data.dlambda_dx_prox) + data.dlambda_dx_prox.size() * sizeof(double);
    total_size += sizeof(data.drhs_prox) + data.drhs_prox.size() * sizeof(double);
    total_size += sizeof(data.jointTorqueRegressor) + data.jointTorqueRegressor.size() * sizeof(double);

    // Matrix6x (6 x Dynamic matrices)
    total_size += sizeof(data.dHdq) + data.dHdq.size() * sizeof(double);
    total_size += sizeof(data.dFdq) + data.dFdq.size() * sizeof(double);
    total_size += sizeof(data.dFdv) + data.dFdv.size() * sizeof(double);
    total_size += sizeof(data.dFda) + data.dFda.size() * sizeof(double);
    total_size += sizeof(data.SDinv) + data.SDinv.size() * sizeof(double);
    total_size += sizeof(data.UDinv) + data.UDinv.size() * sizeof(double);
    total_size += sizeof(data.IS) + data.IS.size() * sizeof(double);
    total_size += sizeof(data.Ag) + data.Ag.size() * sizeof(double);
    total_size += sizeof(data.dAg) + data.dAg.size() * sizeof(double);
    total_size += sizeof(data.J) + data.J.size() * sizeof(double);
    total_size += sizeof(data.dJ) + data.dJ.size() * sizeof(double);
    total_size += sizeof(data.ddJ) + data.ddJ.size() * sizeof(double);
    total_size += sizeof(data.psid) + data.psid.size() * sizeof(double);
    total_size += sizeof(data.psidd) + data.psidd.size() * sizeof(double);
    total_size += sizeof(data.dVdq) + data.dVdq.size() * sizeof(double);
    total_size += sizeof(data.dAdq) + data.dAdq.size() * sizeof(double);
    total_size += sizeof(data.dAdv) + data.dAdv.size() * sizeof(double);

    // RowMatrixXs (Row-major dynamic matrices)
    total_size += sizeof(data.dtau_dq) + data.dtau_dq.size() * sizeof(double);
    total_size += sizeof(data.dtau_dv) + data.dtau_dv.size() * sizeof(double);
    total_size += sizeof(data.ddq_dq) + data.ddq_dq.size() * sizeof(double);
    total_size += sizeof(data.ddq_dv) + data.ddq_dv.size() * sizeof(double);
    total_size += sizeof(data.ddq_dtau) + data.ddq_dtau.size() * sizeof(double);

    // Matrix3x and RowVectorXs
    total_size += sizeof(data.Jcom) + data.Jcom.size() * sizeof(double);
    total_size += sizeof(data.staticRegressor) + data.staticRegressor.size() * sizeof(double);
    total_size += sizeof(data.kineticEnergyRegressor) + data.kineticEnergyRegressor.size() * sizeof(double);
    total_size += sizeof(data.potentialEnergyRegressor) + data.potentialEnergyRegressor.size() * sizeof(double);

    // BodyRegressorType (6x10 matrix)
    total_size += sizeof(data.bodyRegressor);

    // Vector3 vectors (PINOCCHIO_ALIGNED_STD_VECTOR(Vector3))
    total_size += sizeof(data.com) + data.com.capacity() * sizeof(decltype(data.com)::value_type);
    total_size += sizeof(data.vcom) + data.vcom.capacity() * sizeof(decltype(data.vcom)::value_type);
    total_size += sizeof(data.acom) + data.acom.capacity() * sizeof(decltype(data.acom)::value_type);

    // Fixed-size Matrix6 types
    total_size += sizeof(data.Itmp);
    total_size += sizeof(data.M6tmp);
    total_size += sizeof(data.M6tmpR);
    total_size += sizeof(data.M6tmpR2);

    // Fixed-size Force and Inertia types
    total_size += sizeof(data.hg);
    total_size += sizeof(data.dhg);
    total_size += sizeof(data.Ig);

    // Scalar types
    total_size += sizeof(data.kinetic_energy);
    total_size += sizeof(data.potential_energy);
    total_size += sizeof(data.mechanical_energy);

    // std::vector<Scalar> mass
    total_size += sizeof(data.mass) + data.mass.capacity() * sizeof(double);

    // std::vector<int> containers
    total_size += sizeof(data.lastChild) + data.lastChild.capacity() * sizeof(int);
    total_size += sizeof(data.nvSubtree) + data.nvSubtree.capacity() * sizeof(int);
    total_size += sizeof(data.start_idx_v_fromRow) + data.start_idx_v_fromRow.capacity() * sizeof(int);
    total_size += sizeof(data.end_idx_v_fromRow) + data.end_idx_v_fromRow.capacity() * sizeof(int);
    total_size += sizeof(data.idx_vExtended_to_idx_v_fromRow) + data.idx_vExtended_to_idx_v_fromRow.capacity() * sizeof(int);
    total_size += sizeof(data.parents_fromRow) + data.parents_fromRow.capacity() * sizeof(int);
    total_size += sizeof(data.mimic_parents_fromRow) + data.mimic_parents_fromRow.capacity() * sizeof(int);
    total_size += sizeof(data.non_mimic_parents_fromRow) + data.non_mimic_parents_fromRow.capacity() * sizeof(int);
    total_size += sizeof(data.nvSubtree_fromRow) + data.nvSubtree_fromRow.capacity() * sizeof(int);
    total_size += sizeof(data.par_cons_ind) + data.par_cons_ind.capacity() * sizeof(int);
    total_size += sizeof(data.constraint_ind) + data.constraint_ind.capacity() * sizeof(int);
    total_size += sizeof(data.constraints_supported_dim) + data.constraints_supported_dim.capacity() * sizeof(int);

    // std::vector<JointIndex>
    total_size += sizeof(data.mimic_subtree_joint) + data.mimic_subtree_joint.capacity() * sizeof(pinocchio::JointIndex);

    // std::vector<size_t> containers
    total_size += sizeof(data.accumulation_descendant) + data.accumulation_descendant.capacity() * sizeof(size_t);
    total_size += sizeof(data.accumulation_ancestor) + data.accumulation_ancestor.capacity() * sizeof(size_t);
    total_size += sizeof(data.joints_supporting_constraints) + data.joints_supporting_constraints.capacity() * sizeof(size_t);
    total_size += sizeof(data.accumulation_joints) + data.accumulation_joints.capacity() * sizeof(size_t);

    // Nested vector containers
    total_size += sizeof(data.supports_fromRow);
    for (const auto& support_vec : data.supports_fromRow) {
        total_size += sizeof(support_vec) + support_vec.capacity() * sizeof(int);
    }

    total_size += sizeof(data.constraints_on_joint);
    for (const auto& constraint_vec : data.constraints_on_joint) {
        total_size += sizeof(constraint_vec) + constraint_vec.capacity() * sizeof(size_t);
    }

    // std::set containers
    total_size += sizeof(data.constraints_supported);
    for (const auto& constraint_set : data.constraints_supported) {
        total_size += sizeof(constraint_set) + constraint_set.size() * (sizeof(size_t) + sizeof(void*) * 3); // approximate overhead for std::set nodes
    }

    // PINOCCHIO_ALIGNED_STD_VECTOR of complex types
    total_size += sizeof(data.Fcrb) + data.Fcrb.capacity() * sizeof(decltype(data.Fcrb)::value_type);
    total_size += sizeof(data.KA) + data.KA.capacity() * sizeof(decltype(data.KA)::value_type);

    // PINOCCHIO_ALIGNED_STD_VECTOR(MatrixXs) - need to account for individual matrix sizes
    total_size += sizeof(data.LA);
    for (const auto& matrix : data.LA) {
        total_size += sizeof(matrix) + matrix.size() * sizeof(double);
    }

    total_size += sizeof(data.KAS);
    for (const auto& matrix : data.KAS) {
        total_size += sizeof(matrix) + matrix.size() * sizeof(double);
    }

    // PINOCCHIO_ALIGNED_STD_VECTOR(VectorXs)
    total_size += sizeof(data.lA);
    for (const auto& vector : data.lA) {
        total_size += sizeof(vector) + vector.size() * sizeof(double);
    }

    total_size += sizeof(data.lambdaA);
    for (const auto& vector : data.lambdaA) {
        total_size += sizeof(vector) + vector.size() * sizeof(double);
    }

    // Eigen LLT decomposition objects
    total_size += sizeof(data.llt_JMinvJt);
    total_size += sizeof(data.osim_llt);

    // Tensor3x objects (3D tensors)
    total_size += sizeof(data.kinematic_hessians);
    // Note: Tensor size calculation might need more specific handling depending on Tensor implementation

    total_size += sizeof(data.d2tau_dqdq);
    total_size += sizeof(data.d2tau_dvdv);
    total_size += sizeof(data.d2tau_dqdv);
    total_size += sizeof(data.d2tau_dadq);

    // ContactCholeskyDecomposition
    total_size += sizeof(data.contact_chol);

    return total_size;
}

int main()
{
    std::string urdf_path = "/home/quant/research/Galileo/resources/go1/urdf/go1.urdf";

    Model_t model = Model_t();
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<double, 0>(), model);
    Data_t data = Data_t(model);

    auto q = pinocchio::neutral(model);
    auto v = Eigen::VectorXd::Zero(model.nv);
    pinocchio::forwardKinematics(model, data, q, v);
    pinocchio::updateFramePlacements(model, data);
    pinocchio::computeAllTerms(model, data, q, v);

    std::cout << "Model sizeof (stack only): " << sizeof(model) << " bytes" << std::endl;
    std::cout << "Model total memory footprint: " << get_model_size(model) << " bytes" << std::endl;

    std::cout << "\nData sizeof (stack only): " << sizeof(data) << " bytes" << std::endl;
    std::cout << "Data total memory footprint: " << get_data_size(data) << " bytes" << std::endl;

    std::cout << "\nCombined total memory footprint: " << (get_model_size(model) + get_data_size(data)) << " bytes" << std::endl;

    // Print some additional details for verification
    std::cout << "\nModel details:" << std::endl;
    std::cout << "  Number of joints: " << model.njoints << std::endl;
    std::cout << "  Number of bodies: " << model.nbodies << std::endl;
    std::cout << "  Number of frames: " << model.nframes << std::endl;
    std::cout << "  Configuration dimension (nq): " << model.nq << std::endl;
    std::cout << "  Velocity dimension (nv): " << model.nv << std::endl;

    return 0;
}
