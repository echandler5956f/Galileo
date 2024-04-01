#include "atlas_test.h"

int main(int argc, char **argv)
{
    std::cout << "Starting Atlas Test" << std::endl;
    Eigen::VectorXd q0_vec = Eigen::VectorXd::Zero(37);
    q0_vec(2) = 0.864977;
    q0_vec(6) = 1.;

    q0_vec(11) = -1.5;
    q0_vec(27) = -0.5;
    q0_vec(28) = 1.;
    q0_vec(29) = -0.5;

    q0_vec(19) = 1.5;
    q0_vec(33) = -0.5;
    q0_vec(34) = 1.;
    q0_vec(35) = -0.5;

    std::cout << "q0_vec: " << q0_vec << std::endl;

    // std::vector<int> knot_num = {25};
    // std::vector<double> knot_time = {0.25};
    // std::vector<uint> mask_vec = {0b11}; // standing

    // std::vector<int> knot_num = {30, 30, 30};
    // std::vector<double> knot_time = {0.05, 1.0, 0.05};
    // std::vector<uint> mask_vec = {0b11, 0b00, 0b11}; // static walk

    std::vector<int> knot_num = {25, 16, 16, 16, 16, 25};
    std::vector<double> knot_time = {0.25, 0.1667, 0.1667, 0.1667, 0.1667, 0.25};
    std::vector<uint> mask_vec = {0b11, 0b10, 0b01, 0b10, 0b01, 0b11}; // static walk

    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;
    std::vector<galileo::legged::environment::SurfaceID> tmp_surfaces;

    for (int i = 0; i < mask_vec.size(); i++)
    {
        tmp_surfaces.clear();
        for (int j = 0; j < end_effector_names.size(); j++)
        {
            if (mask_vec[i] & (1 << j))
            {
                tmp_surfaces.push_back(0);
                std::cout << "Adding surface 0" << std::endl;
            }
            else
            {
                tmp_surfaces.push_back(galileo::legged::environment::NO_SURFACE);
                std::cout << "Adding surface -1" << std::endl;
            }
        }
        contact_surfaces.push_back(tmp_surfaces);
    }

    galileo::legged::LeggedInterface solver_interface;

    solver_interface.LoadModel(atlas_location, end_effector_names);
    solver_interface.LoadParameters(atlas_parameter_location);

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
    Xf(solver_interface.states()->q_index) = 0.4;

    solver_interface.Initialize(X0, Xf);
    solver_interface.Update(X0, Xf);

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., solver_interface.getRobotModel()->contact_sequence->getDT());
    Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
    Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());
    // solver_interface.GetSolution(new_times, new_states, new_inputs);
    solver_interface.VisualizeSolutionAndConstraints(new_times, new_states, new_inputs);

    return 0;
}