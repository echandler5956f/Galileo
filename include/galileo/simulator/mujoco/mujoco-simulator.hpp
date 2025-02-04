#ifndef __galileo_simulator_mujoco_mujoco_simulator_hpp__
#define __galileo_simulator_mujoco_mujoco_simulator_hpp__

#include "galileo/simulator/simulator-base.hpp"

#include <mujoco/mujoco.h>
#include <GLFW/glfw3.h>
#include <memory>
#include <cstring>

namespace galileo
{
    namespace simulator
    {

        namespace mujoco
        {

            static GLFWwindow *window;
            static mjvCamera cam;  // abstract camera
            static mjvOption opt;  // visualization options
            static mjvScene scn;   // abstract scene
            static mjrContext con; // custom GPU context

            // moue interaction
            static bool button_left = false;
            static bool button_middle = false;
            static bool button_right = false;
            static double lastx = 0;
            static double lasty = 0;

            static mjModel *m; // MuJoCo model
            static mjData *d;  // MuJoCo data

            static mjfGeneric user_loop_func; // user-defined loop function

            // keyboard callback
            void keyboard(GLFWwindow *window, int key, int scancode, int act, int mods)
            {
                // backspace: reset simulation
                if (act == GLFW_PRESS && key == GLFW_KEY_BACKSPACE)
                {
                    mj_resetData(m, d);
                    mj_forward(m, d);
                }
            }

            // mouse button callback
            void mouse_button(GLFWwindow *window, int button, int act, int mods)
            {
                // update button state
                button_left = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) ==
                               GLFW_PRESS);
                button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) ==
                                 GLFW_PRESS);
                button_right = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) ==
                                GLFW_PRESS);

                // update mouse position
                glfwGetCursorPos(window, &lastx, &lasty);
            }

            // mouse move callback
            void mouse_move(GLFWwindow *window, double xpos, double ypos)
            {
                // no buttons down: nothing to do
                if (!button_left && !button_middle && !button_right)
                {
                    return;
                }

                // compute mouse displacement, save
                double dx = xpos - lastx;
                double dy = ypos - lasty;
                lastx = xpos;
                lasty = ypos;

                // get current window size
                int width, height;
                glfwGetWindowSize(window, &width, &height);

                // get shift key state
                bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                                  glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

                // determine action based on mouse button
                mjtMouse action;
                if (button_right)
                {
                    action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
                }
                else if (button_left)
                {
                    action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
                }
                else
                {
                    action = mjMOUSE_ZOOM;
                }

                // move camera
                mjv_moveCamera(m, action, dx / height, dy / height, &scn, &cam);
            }

            // scroll callback
            void scroll(GLFWwindow *window, double xoffset, double yoffset)
            {
                // emulate vertical mouse motion = 5% of window height
                mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05 * yoffset, &scn, &cam);
            }

            void dummy_loop() {}

            /**
             * Must be called after MujocoSimulator().
             */
            void mjcontroller(const mjModel *m, mjData *d)
            {
                user_loop_func(m, d);
            }

            class MujocoSimulator : public SimulatorBase
            {
            public:
                MujocoSimulator(mjfGeneric loop_func, std::function<void(void)> exit_func = []() {}) : SimulatorBase(std::bind(dummy_loop), exit_func)
                {
                    user_loop_func = loop_func;
                }

                MujocoSimulator(const MujocoSimulator &) = delete;

                MujocoSimulator &operator=(const MujocoSimulator &) = delete;

                ~MujocoSimulator() = default;

                void Initialize(const std::string &model_path) override
                {
                    char error[1000] = "Could not load XML model";
                    m = mj_loadXML(model_path.c_str(), nullptr, error, 1000);
                    if (!m)
                    {
                        std::cerr << "Load model error: " << error << std::endl;
                        std::exit(1);
                    }
                    d = mj_makeData(m);
                    initializeImp();
                }

                void Initialize(mjModel *model, mjData *data)
                {
                    m = model;
                    d = data;
                    initializeImp();
                }

                void Configure(void *config) override
                {
                }

                void Finalize() override
                {
                }

                void LoopPrologue() override
                {
                    mjtNum simstart = d->time;
                    while (d->time - simstart < 1.0 / 60.0)
                    {
                        mj_step(m, d);
                    }

                    mjrRect viewport = {0, 0, 0, 0};
                    glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

                    // update scene and render
                    mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
                    mjr_render(viewport, &scn, &con);

                    // swap OpenGL buffers (blocking call due to v-sync)
                    glfwSwapBuffers(window);

                    // process pending GUI events, call GLFW callbacks
                    glfwPollEvents();
                }

                void LoopEpilogue() override
                {
                }

                void SetMotorValue(int id, double u)
                {
                    d->ctrl[id] = u;
                }

                double GetTime() const
                {
                    return d->time;
                }

                void Exit() override
                {
                    mjv_freeScene(&scn);
                    mjr_freeContext(&con);

                    // free MuJoCo model and data
                    mj_deleteData(d);
                    mj_deleteModel(m);

                    // terminate GLFW (crashes with Linux NVidia drivers)
#if defined(__APPLE__) || defined(_WIN32)
                    glfwTerminate();
#endif
                }

            public:
                void initializeImp()
                {
                    // init GLFW
                    if (!glfwInit())
                    {
                        mju_error("Could not initialize GLFW");
                    }

                    // create window, make OpenGL context current, request v-sync
                    window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
                    glfwMakeContextCurrent(window);
                    glfwSwapInterval(1);

                    // initialize visualization data structures
                    mjv_defaultCamera(&cam);
                    mjv_defaultOption(&opt);
                    mjv_defaultScene(&scn);
                    mjr_defaultContext(&con);

                    // create scene and context
                    mjv_makeScene(m, &scn, 2000);
                    mjr_makeContext(m, &con, mjFONTSCALE_150);

                    // install GLFW mouse and keyboard callbacks
                    glfwSetKeyCallback(window, keyboard);
                    glfwSetCursorPosCallback(window, mouse_move);
                    glfwSetMouseButtonCallback(window, mouse_button);
                    glfwSetScrollCallback(window, scroll);
                    mjcb_control = mjcontroller;
                }

            }; // class MujocoSimulator

        } // namespace mujoco

    } // namespace simulator

} // namespace galileo

#endif // __galileo_simulator_mujoco_mujoco_simulator_hpp__