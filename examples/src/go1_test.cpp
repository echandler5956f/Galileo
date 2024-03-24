#include "go1_test.h"

int main(int argc, char **argv)
{
    ConfigVector q0_vec = (ConfigVector(19) << 0., 0., 0.339, 0., 0., 0., 1., 0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3).finished();

    // std::vector<int> knot_num = {5, 30, 30, 5};
    // std::vector<double> knot_time = {0.05, 0.3, 0.3, 0.05};
    // std::vector<uint> mask_vec = {0b1111, 0b1001, 0b0110,  0b1111}; // trot

    std::vector<int> knot_num = {5, 30, 30, 30, 30, 5};
    std::vector<double> knot_time = {0.05, 0.3, 0.3, 0.3, 0.3, 0.05};
    std::vector<uint> mask_vec = {0b1111, 0b1001, 0b0110, 0b1001, 0b0110, 0b1111}; // trot

    // std::vector<int> knot_num = {5, 30};
    // std::vector<double> knot_time = {0.05, 0.3};
    // std::vector<uint> mask_vec = {0b1111, 0b1001}; // trot

    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;
    std::vector<galileo::legged::environment::SurfaceID> tmp_surfaces;

    for (int i = 0; i < mask_vec.size(); i++)
    {
        tmp_surfaces.clear();
        for (int j = 0; j < 4; j++)
        {
            if (mask_vec[i] & (8 >> j))
            {
                tmp_surfaces.push_back(0);
            }
            else
            {
                tmp_surfaces.push_back(galileo::legged::environment::NO_SURFACE);
            }
        }
        contact_surfaces.push_back(tmp_surfaces);
    }

    galileo::legged::LeggedInterface solver_interface;

    solver_interface.LoadModel(go1_location, end_effector_names);
    solver_interface.LoadParameters(go1_parameter_location);

    casadi::DM X0 = casadi::DM::zeros(solver_interface.states()->nx, 1);
    int q0_idx = solver_interface.states()->q_index;
    for (int j = 0; j < solver_interface.states()->nq; j++)
    {
        X0(q0_idx + j) = q0_vec(j);
    }

    solver_interface.addSurface(environment::createInfiniteGround());

    solver_interface.setContactSequence(
        knot_num,
        knot_time,
        mask_vec,
        contact_surfaces);

    casadi::DM Xf = X0;
    Xf(solver_interface.states()->q_index) = 0.2;
    // Xf(solver_interface.states()->q_index) = 0.2 * sqrt(2);
    // Xf(solver_interface.states()->q_index + 2) = 0.2 * sqrt(2);
    // Xf(solver_interface.states()->q_index + 3) = 0.;
    // Xf(solver_interface.states()->q_index + 4) = 0.;
    // Xf(solver_interface.states()->q_index + 5) = 0.3826834;
    // Xf(solver_interface.states()->q_index + 6) = 0.9238795;

    solver_interface.Initialize(X0, Xf);
    solver_interface.Update(X0, Xf);

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., solver_interface.getRobotModel()->contact_sequence->getDT());
    Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
    Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());
    // solver_interface.GetSolution(new_times, new_states, new_inputs);
    solver_interface.VisualizeSolutionAndConstraints(new_times, new_states, new_inputs);

    return 0;
}