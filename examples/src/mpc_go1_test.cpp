#include "go1_test.h"

int main(int argc, char **argv)
{
    std::vector<std::string> end_effector_names;
    std::vector<int> knot_num;
    std::vector<double> knot_time;
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;

    std::vector<double> q0_vec;
    std::vector<double> qf_vec;

    galileo::legged::helper::ReadProblemFromParameterFile(problem_parameter_location,
                                                          end_effector_names,
                                                          knot_num,
                                                          knot_time,
                                                          contact_surfaces,
                                                          q0_vec,
                                                          qf_vec);

    galileo::legged::LeggedInterface solver_interface;

    solver_interface.LoadModel(robot_location, end_effector_names);
    solver_interface.LoadParameters(solver_parameter_location);

    int nx = solver_interface.states()->nx;
    int q_idx = solver_interface.states()->q_index;

    casadi::DM X0;
    std::vector<double> X0_vec = galileo::legged::helper::getXfromq(solver_interface.states()->nx, q_idx, q0_vec);
    galileo::tools::vectorToCasadi(X0_vec, nx, 1, X0);

    casadi::DM Xf;
    std::vector<double> Xf_vec = galileo::legged::helper::getXfromq(solver_interface.states()->nx, q_idx, qf_vec);
    galileo::tools::vectorToCasadi(Xf_vec, nx, 1, Xf);

    solver_interface.addSurface(environment::createInfiniteGround());

    std::vector<uint> mask_vec = galileo::legged::helper::getMaskVectorFromContactSurfaces(contact_surfaces);

    solver_interface.setContactSequence(
        knot_num,
        knot_time,
        mask_vec,
        contact_surfaces);

    solver_interface.Initialize(X0, Xf, 1, 1.3);

    int mpc_iterations = 100;
    double dt = solver_interface.getSolutionDT() / (mpc_iterations);

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(mpc_iterations, 0., dt * (mpc_iterations));
    Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
    Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());

    auto time_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < mpc_iterations; i++)
    {
        solver_interface.Update(new_times[i], X0, Xf);
        Eigen::VectorXd time = Eigen::VectorXd::Zero(1);
        time(0) = new_times[i] + dt;
        Eigen::MatrixXd state = Eigen::MatrixXd::Zero(solver_interface.states()->nx, 1);
        Eigen::MatrixXd input = Eigen::MatrixXd::Zero(solver_interface.states()->nu, 1);
        solver_interface.GetSolution(time, state, input);
        new_states.col(i) = state;
        new_inputs.col(i) = input;
        tools::eigenToCasadi(state, X0);
    }
    auto time_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = time_end - time_start;
    auto mpc_frequency = mpc_iterations / elapsed_seconds.count();
    std::cout << "MPC frequency: " << mpc_frequency << " Hz" << std::endl;
    solver_interface.PlotTrajectories(new_times, new_states, new_inputs);
    solver_interface.PlotConstraints(new_times, new_states, new_inputs);

    // Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., solver_interface.getRobotModel()->contact_sequence->getDT());
    // Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
    // Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());
    // // solver_interface.GetSolution(new_times, new_states, new_inputs);
    // solver_interface.VisualizeSolutionAndConstraints(new_times, new_states, new_inputs);

    return 0;
}