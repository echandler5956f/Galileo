#include <cassert>
#include <chrono>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Dense>
#include <GLFW/glfw3.h>
#include <mujoco/mujoco.h>

#include "galileo/common/utils/resource_finder.hpp"

// =================================================================================================
// 1. GLOBAL VARIABLES AND DATA STRUCTURES
// =================================================================================================

// --- MuJoCo Globals ---
mjModel *m = nullptr; // MuJoCo model
mjData *d = nullptr;  // MuJoCo data
mjvCamera cam;        // Abstract camera
mjvOption opt;        // Visualization options
mjvScene scn;         // Abstract scene
mjrContext con;       // GL rendering context

// --- GLFW and UI Globals ---
GLFWwindow *window = nullptr;
std::mutex mtx; // Mutex for thread-safe access to MuJoCo data
bool button_left = false;
bool button_middle = false;
bool button_right = false;
double lastx = 0;
double lasty = 0;

// --- Controller Configuration ---
const double CONTROL_FREQUENCY = 200.0;     // Hz
const double SIMULATION_FREQUENCY = 1000.0; // Hz (must be >= control frequency)
const int CONTROL_DECIMATION = SIMULATION_FREQUENCY / CONTROL_FREQUENCY;
const double RENDER_FPS = 60.0;

// PD Gains
double kp = 30.0;
double kd = 1.0;
Eigen::VectorXd q_des; // Desired joint positions

/**
 * @struct RobotState
 * @brief A clean container for the robot's full state using Eigen vectors.
 *
 */
struct RobotState
{
    // Using Dynamic-sized vectors for generality
    Eigen::VectorXd q; // Generalized positions: [base_pos (3), base_quat (4), joint_pos (nq-7)]
    Eigen::VectorXd v; // Generalized velocities: [base_linear_vel (3), base_angular_vel (3), joint_vel (nv-6)]
};

// =================================================================================================
// 2. CORE FUNCTIONS: STATE EXTRACTION & CONTROL
// =================================================================================================

/**
 * @brief Extracts the full state of the robot from MuJoCo's mjData.
 *
 * @param model A const pointer to the mjModel.
 * @param data A const pointer to the mjData.
 * @return A RobotState object containing the current state.
 */
RobotState get_robot_state(const mjModel *model, const mjData *data)
{
    assert(model != nullptr && data != nullptr && "MuJoCo model and data must be initialized.");
    assert(model->nq >= 7 && model->nv >= 6 && "Model must have a floating base (at least 7 qpos and 6 qvel).");

    RobotState state;
    state.q = Eigen::Map<const Eigen::VectorXd>(data->qpos, model->nq);
    state.v = Eigen::Map<const Eigen::VectorXd>(data->qvel, model->nv);

    return state;
}

/**
 * @brief Applies a vector of joint torques to the MuJoCo simulation.
 *
 * @param model A const pointer to the mjModel.
 * @param data A pointer to the mjData where control inputs will be written.
 * @param torques An Eigen::VectorXd of joint torques. Its size must match the number of actuators.
 */
void set_joint_torques(const mjModel *model, mjData *data, const Eigen::VectorXd &torques)
{
    assert(model != nullptr && data != nullptr && "MuJoCo model and data must be initialized.");
    assert(torques.size() == model->nu && "Torque vector size must match the number of actuators.");

    // Eigen::Map allows writing directly from an Eigen object to a C-style array.
    Eigen::Map<Eigen::VectorXd>(data->ctrl, model->nu) = torques;
}

/**
 * @brief Computes joint torques using a Proportional-Derivative (PD) controller.
 *
 * @param state The current RobotState of the system.
 * @param q_desired The desired joint positions (Eigen::VectorXd).
 * @param Kp The proportional gain.
 * @param Kd The derivative gain.
 * @param nq_joints The number of joints (nq - 7 for a floating base).
 * @param nv_joints The number of joint velocities (nv - 6 for a floating base).
 * @return An Eigen::VectorXd containing the computed control torques.
 */
Eigen::VectorXd pd_controller(const RobotState &state, const Eigen::VectorXd &q_desired,
                              double Kp, double Kd, int nq_joints, int nv_joints)
{
    assert(q_desired.size() == nq_joints && "Desired position vector size must match number of joints.");
    assert(state.q.size() >= 7 + nq_joints && "State vector q is too small.");
    assert(state.v.size() >= 6 + nv_joints && "State vector v is too small.");

    Eigen::VectorXd q_joints_current = state.q.segment(7, nq_joints);
    Eigen::VectorXd v_joints_current = state.v.segment(6, nv_joints);

    // Compute position and velocity errors.
    Eigen::VectorXd error_q = q_desired - q_joints_current;
    Eigen::VectorXd error_v = -v_joints_current; // Target velocity is zero.

    // Calculate PD control law.
    Eigen::VectorXd torques = Kp * error_q + Kd * error_v;

    return torques;
}

// =================================================================================================
// 3. GLFW CALLBACKS FOR INTERACTIVE RENDERING
// =================================================================================================

void keyboard(GLFWwindow *window, int key, int scancode, int act, int mods)
{
    if (act == GLFW_PRESS && key == GLFW_KEY_BACKSPACE)
    {
        mj_resetData(m, d);
        mj_forward(m, d);
        std::cout << "Simulation reset." << std::endl;
    }
}

void mouse_button(GLFWwindow *window, int button, int act, int mods)
{
    button_left = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS);
    button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS);
    button_right = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS);
    glfwGetCursorPos(window, &lastx, &lasty);
}

void mouse_move(GLFWwindow *window, double xpos, double ypos)
{
    if (!button_left && !button_middle && !button_right)
    {
        return;
    }

    double dx = xpos - lastx;
    double dy = ypos - lasty;
    lastx = xpos;
    lasty = ypos;

    int width, height;
    glfwGetWindowSize(window, &width, &height);

    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

    if (button_right)
    {
        mod_shift ? mjv_moveCamera(m, mjMOUSE_ZOOM, dx, dy, &scn, &cam)
                  : mjv_moveCamera(m, mjMOUSE_ROTATE_V, dx, dy, &scn, &cam);
    }
    else if (button_left)
    {
        mjv_moveCamera(m, mjMOUSE_ROTATE_H, dx, dy, &scn, &cam);
    }
    else if (button_middle)
    {
        mjv_moveCamera(m, mjMOUSE_MOVE_V, dx, dy, &scn, &cam);
    }
}

void scroll(GLFWwindow *window, double xoffset, double yoffset)
{
    mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -5 * yoffset, &scn, &cam);
}

// =================================================================================================
// 4. MAIN SIMULATION FUNCTION
// =================================================================================================

int main(int argc, char **argv)
{
    // --------------------------- 4.1. INITIALIZATION ---------------------------
    std::cout << "Initializing MuJoCo PD Controller..." << std::endl;

    // --- Provide path to your model here ---
    const std::string model_path = galileo::utils::get_resource_path("go2/mjx_scene_force.xml");

    // Load the model
    char error[1000] = "Could not load model file";
    m = mj_loadXML(model_path.c_str(), nullptr, error, 1000);
    assert(m && "Failed to load model. Check the path and XML syntax.");
    if (!m)
    {
        std::cerr << "Error loading model: " << error << std::endl;
        return -1;
    }

    // Make data
    d = mj_makeData(m);
    assert(d && "Failed to create MuJoCo data.");

    int home_id = mj_name2id(m, mjOBJ_KEY, "home");

    for (size_t i = 0; i < m->nq; i++)
    {
        d->qpos[i] = m->key_qpos[home_id * m->nq + i];
    }

    // Initialize desired joint positions to the model's home configuration
    Eigen::VectorXd q_home = Eigen::VectorXd::Zero(m->nq);
    for (size_t i = 0; i < m->nq; i++)
    {
        q_home(i) = d->qpos[i];
    }

    int nq_joints = m->nq - 7;
    q_des = q_home.segment(7, nq_joints);

    std::cout << "Model loaded successfully. NQ: " << m->nq << ". NV: " << m->nv << ". NU: " << m->nu << "." << std::endl;
    std::cout << "Initial desired joint positions set to model home." << std::endl;

    std::cout << "Home: " << q_home.transpose() << std::endl;
    std::cout << "Desired: " << q_des.transpose() << std::endl;

    // Initialize GLFW
    if (!glfwInit())
    {
        std::cerr << "Could not initialize GLFW." << std::endl;
        return -1;
    }
    window = glfwCreateWindow(1200, 900, "MuJoCo PD Control Demo", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize MuJoCo visualization
    mjv_defaultCamera(&cam);
    mjv_defaultOption(&opt);
    mjv_defaultScene(&scn);
    mjr_defaultContext(&con);
    mjv_makeScene(m, &scn, 2000);
    mjr_makeContext(m, &con, mjFONTSCALE_150);

    // Install GLFW callbacks
    glfwSetKeyCallback(window, keyboard);
    glfwSetCursorPosCallback(window, mouse_move);
    glfwSetMouseButtonCallback(window, mouse_button);
    glfwSetScrollCallback(window, scroll);

    // --------------------------- 4.2. MAIN SIMULATION LOOP ---------------------------
    std::cout << "Starting simulation loop..." << std::endl;

    auto last_control_time = std::chrono::high_resolution_clock::now();
    auto last_render_time = last_control_time;
    const auto control_dt = std::chrono::duration<double>(1.0 / CONTROL_FREQUENCY);
    const auto render_dt = std::chrono::duration<double>(1.0 / RENDER_FPS);

    auto current_time = std::chrono::high_resolution_clock::now();

    while (!glfwWindowShouldClose(window))
    {
        current_time = std::chrono::high_resolution_clock::now();

        // --- Control Step (runs at lower frequency) ---
        if ((current_time - last_control_time) >= control_dt)
        {
            // Use a mutex to ensure data is not being written by mj_step while we read it
            std::lock_guard<std::mutex> lock(mtx);

            // Step 2: Extract current robot state
            RobotState current_state = get_robot_state(m, d);

            // Step 3: Compute and apply control torques
            int nv_joints = m->nv - 6;
            Eigen::VectorXd torques = pd_controller(current_state, q_des, kp, kd, nq_joints, nv_joints);
            set_joint_torques(m, d, torques);

            last_control_time = std::chrono::high_resolution_clock::now();
        }

        // --- Physics Step (runs at higher frequency) ---
        {
            std::lock_guard<std::mutex> lock(mtx);
            mj_step(m, d);
        }

        // --- Rendering Step ---
        mjrRect viewport = {0, 0, 0, 0};
        current_time = std::chrono::high_resolution_clock::now();
        if ((current_time - last_render_time) >= render_dt)
        {
            glfwGetFramebufferSize(window, &viewport.width, &viewport.height);
            mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
            mjr_render(viewport, &scn, &con);
            last_render_time = std::chrono::high_resolution_clock::now();
        }

        // Display simulation info
        char info_str[1024];
        sprintf(info_str, "Time: %.2f s\nFrequency: %.1f Hz", d->time, 1.0 / (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - last_control_time).count()));
        mjr_text(mjFONT_NORMAL, info_str, &con, 10, viewport.height - 20, 0, 0, 0);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // --------------------------- 4.3. CLEANUP ---------------------------
    std::cout << "Simulation finished. Cleaning up resources." << std::endl;
    mjv_freeScene(&scn);
    mjr_freeContext(&con);
    mj_deleteData(d);
    mj_deleteModel(m);
    glfwTerminate();

    return 0;
}
