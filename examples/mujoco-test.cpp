#include <iostream>
#include <galileo/simulator/mujoco/mujoco-simulator.hpp>
// #include <galileo/reactive/wbc-base4.hpp>
#include <galileo/reactive/wbc-base.hpp>
#include <pinocchio/algorithm/rnea.hpp>

#include <chrono>
#include <thread>

using namespace galileo;

typedef pinocchio::ModelTpl<double> ModelType;
typedef pinocchio::DataTpl<double> DataType;

void loop(const mjModel *m, mjData *d);

simulator::mujoco::MujocoSimulator mjEnv(loop);
reactive::WBCBase<double> wbc_base;
// reactive::WBCBase wbc_base;

ModelType model;
DataType data;

std::vector<EndEffector> ees;
Eigen::Vector<double, 12> static_force;

Eigen::VectorXd state_desired;
Eigen::VectorXd control_desired;

bool controller_started = false;

int main(int argc, char *argv[])
{
    std::string model_path = "/home/quant/dial_mpc_ws/src/dial-mpc/models/unitree_go2/mjx_scene_force.xml";
    // std::string pin_model_path = "/home/quant/dial_mpc_ws/src/dial-mpc/models/unitree_go2/mjx_go2_force.xml";
    std::string pin_model_path = "/home/quant/dial_mpc_ws/src/dial-mpc/models/unitree_go2/go2_description.urdf";

    model = ModelType();
    // pinocchio::mjcf::buildModel(pin_model_path, pinocchio::JointModelFreeFlyer(), model);
    pinocchio::urdf::buildModel(pin_model_path, pinocchio::JointModelFreeFlyer(), model);
    data = DataType(model);

    for (int i = 0; i < model.njoints; i++)
    {
        std::cout << "pin joint name: " << model.names[i] << std::endl;
    }

    std::cout << "nq: " << model.nq << std::endl;
    std::cout << "nv: " << model.nv << std::endl;

    std::vector<std::string> ee_names = {"FR_foot", "FL_foot", "RR_foot", "RL_foot"};

    reactive::Info info;
    info.nq = model.nq;
    info.nv = model.nv;
    info.actuatedDofNum = 12;
    info.numThreeDofContacts = ee_names.size();

    // std::vector<EndEffector> ees;
    for (int i = 0; i < info.numThreeDofContacts; i++)
    {
        EndEffector ee;
        ee.frame_name = ee_names[i];
        ee.frame_idx = model.getFrameId(ee.frame_name);
        ees.push_back(ee);
    }

    Eigen::VectorXd q_des(info.nq);
    q_des << 0, 0, 0.27, 0, 0, 0, 1, 0, 0.9, -1.8, 0, 0.9, -1.8, 0, 0.9, -1.8, 0, 0.9, -1.8;
    // q_des << 0, 0, 0.445, 0, 0, 0, 1, 0, 0.9, -1.8, 0, 0.9, -1.8, 0, 0.9, -1.8, 0, 0.9, -1.8;

    Eigen::VectorXd v_des = Eigen::VectorXd::Zero(info.nv);
    state_desired = Eigen::VectorXd::Zero(info.nq + info.nv);
    state_desired.head(info.nq) = q_des;
    state_desired.tail(info.nv) = v_des;

    pinocchio::computeTotalMass(model, data);
    pinocchio::forwardKinematics(model, data, q_des, v_des);
    pinocchio::updateFramePlacements(model, data);
    pinocchio::computeJointJacobians(model, data, q_des);

    Eigen::VectorXd static_force_per_leg(3);
    static_force_per_leg << 0, 0, 9.81 * data.mass[0] / 4;

    static_force = Eigen::Vector<double, 12>::Zero();
    for (size_t i = 0; i < info.numThreeDofContacts; ++i)
    {
        static_force.segment<3>(3 * i) = static_force_per_leg;
    }

    // use Jacobian to get the contribution of the static forces to the torques
    Eigen::Matrix<double, 12, 12> J = Eigen::Matrix<double, 12, 12>::Zero();
    for (size_t i = 0; i < info.numThreeDofContacts; ++i)
    {
        Eigen::Matrix<double, 6, 18> jac = Eigen::Matrix<double, 6, 18>::Zero(6, info.nv);
        pinocchio::getFrameJacobian(model, data, ees[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED, jac);
        Eigen::Matrix<double, 3, 18> jac_pos = jac.topRows(3);
        J.block(3 * i, 0, 3, info.actuatedDofNum) = jac_pos.rightCols(info.actuatedDofNum);
    }

    Eigen::Vector<double, 12> tau_des = -J.transpose() * static_force;
    std::cout << "tau desired: " << tau_des.transpose() << std::endl;

    // control is u = [a, f, tau]
    control_desired = Eigen::VectorXd::Zero(info.nv + 3 * info.numThreeDofContacts + info.actuatedDofNum);
    control_desired.head(info.nv) = Eigen::VectorXd::Zero(info.nv);
    control_desired.segment(info.nv, 3 * info.numThreeDofContacts) = static_force;
    control_desired.tail(info.actuatedDofNum) = tau_des;

    std::cout << "control desired: " << control_desired.transpose() << std::endl;

    std::vector<std::string> joint_names = {"FR_hip_joint", "FR_thigh_joint", "FR_calf_joint", "FL_hip_joint", "FL_thigh_joint", "FL_calf_joint", "RR_hip_joint", "RR_thigh_joint", "RR_calf_joint", "RL_hip_joint", "RL_thigh_joint", "RL_calf_joint"};

    // wbc_base = reactive::WBCBase<double>(model, info, ees);
    wbc_base = reactive::WBCBase(model, info, ees);

    std::cout << "Created controller" << std::endl;

    mjEnv.Initialize(model_path);

    // for (int i = 0; i < joint_names.size(); i++)
    // {
    //     std::cout << "mj joint name: " << mj_name2id(simulator::mujoco::m, mjOBJ_JOINT, joint_names[i].c_str()) << std::endl;
    // }

    std::cout << "Starting controller" << std::endl;

    int home_id = mj_name2id(simulator::mujoco::m, mjOBJ_KEY, "home");
    if (home_id < 0)
    {
        std::cerr << "Keyframe 'home' not found" << std::endl;
    }
    mju_copy(simulator::mujoco::d->qpos, simulator::mujoco::m->key_qpos + simulator::mujoco::m->nq * home_id, simulator::mujoco::m->nq);

    mjEnv.Finalize();

    controller_started = true;

    mjEnv.Loop();

    mjEnv.Exit();

    return 0;
}

void loop(const mjModel *m, mjData *d)
{
    if (!controller_started)
    {
        return;
    }

    // Get the measured state
    Eigen::VectorXd q = Eigen::VectorXd::Map(d->qpos, m->nq);
    Eigen::VectorXd v = Eigen::VectorXd::Map(d->qvel, m->nv);

    Eigen::VectorXd state = Eigen::VectorXd::Zero(m->nq + m->nv);
    state.head(m->nq) = q;
    state.tail(m->nv) = v;

    // we need to swap the quaternion positions, since pinocchio uses xyzw and mujoco uses wxyz
    // state from d->qpos is ordered as [posx, posy, posz, quatw, quatx, quaty, quatz, ...]
    // but we need it to be ordered as [posx, posy, posz, quatx, quaty, quatz, quatw, ...] for pinocchio

    Eigen::VectorXd state_swapped = state;
    // Extract the quaternion (of size 4) starting at index 3.
    Eigen::Vector4d quat = state_swapped.segment(3, 4);
    // Create a new quaternion with the desired ordering.
    // That is, if the original is [w, x, y, z], we want [x, y, z, w].
    Eigen::Vector4d quat_swapped;
    quat_swapped << quat(1), quat(2), quat(3), quat(0);
    // Replace the segment in state_swapped.
    state_swapped.segment(3, 4) = quat_swapped;

    // simple PD controller:
    Eigen::VectorXd tau_fb = Eigen::VectorXd::Zero(m->nu);
    Eigen::VectorXd q_des = state_desired.segment(7, 12);
    Eigen::VectorXd v_des = state_desired.tail(12);
    Eigen::VectorXd q_err = q_des - q.tail(12);
    Eigen::VectorXd v_err = v_des - v.tail(12);
    double kp = 50.0;
    double kd = 0.0;
    for (int i = 0; i < m->nu; i++)
    {
        tau_fb[i] = kp * q_err[i] + kd * v_err[i];
    }

    // // tau = control_desired.tail(m->nu);

    pinocchio::forwardKinematics(model, data, state_swapped.head(19), v);
    pinocchio::updateFramePlacements(model, data);
    pinocchio::computeJointJacobians(model, data, state_swapped.head(19));

    Eigen::Matrix<double, 12, 12> J = Eigen::Matrix<double, 12, 12>::Zero();
    for (size_t i = 0; i < 4; ++i)
    {
        Eigen::Matrix<double, 6, 18> jac = Eigen::Matrix<double, 6, 18>::Zero(6, 18);
        pinocchio::getFrameJacobian(model, data, ees[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED, jac);
        Eigen::Matrix<double, 3, 18> jac_pos = jac.topRows(3);
        J.block(3 * i, 0, 3, 12) = jac_pos.rightCols(12);
    }

    Eigen::Vector<double, 12> tau_ff = -J.transpose() * static_force;

    control_desired.tail(m->nu) = tau_ff;

    // control_desired = Eigen::VectorXd::Zero(m->nu + 3 * 4 + 12);
    
    // auto start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXd control_new = wbc_base.update(state_desired, control_desired, state_swapped, 0, 0.0);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;

    std::cout << "state: " << q.transpose() << std::endl;
    std::cout << "control_desired: " << control_desired.transpose() << std::endl;
    std::cout << "control_new: " << control_new.transpose() << std::endl;

    // Eigen::VectorXd tau_new = control_new.tail(m->nu);
    Eigen::VectorXd tau_new = control_new.tail(m->nu) + tau_fb;

    // Set the control
    for (int i = 0; i < m->nu; i++)
    {
        d->ctrl[i] = tau_new[i];
    }

    // std::this_thread::sleep_for(std::chrono::milliseconds(50));
}