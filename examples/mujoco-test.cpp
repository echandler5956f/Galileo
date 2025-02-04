#include <iostream>
#include <galileo/simulator/mujoco/mujoco-simulator.hpp>

using namespace galileo;

void loop(const mjModel *m, mjData *d);

simulator::mujoco::MujocoSimulator mjEnv(loop);

bool controller_started = false;
double mj_start_time = 0.0;
double curr_time = 0.0;

int main(int argc, char *argv[])
{
    std::string model_path = "/home/quant/dial_mpc_ws/src/dial-mpc/models/unitree_go2/mjx_scene_force.xml";

    mjEnv.Initialize(model_path);
    controller_started = true;
    mj_start_time = mjEnv.GetTime();
    mjEnv.Finalize();

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

    // // Set the control
    // for (int i = 0; i < m->nu; i++)
    // {
    //     d->ctrl[i] = tau[i];
    // }
}