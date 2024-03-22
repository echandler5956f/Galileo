#include "go1_test.h"

int main(int argc, char **argv)
{
    ConfigVector q0_vec = (ConfigVector(19) << 0., 0., 0.339, 0., 0., 0., 1., 0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3).finished();

    std::vector<int> knot_num = {30, 30, 30, 30, 30, 30};
    std::vector<double> knot_time = {0.05, 0.3, 0.3, 0.3, 0.3, 0.05};
    std::vector<uint> mask_vec = {0b1111, 0b1001, 0b0110, 0b1001, 0b0110, 0b1111}; // trot

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

    solver_interface.Initialize(X0, X0);
    solver_interface.Update(X0, X0);

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., solver_interface.getRobotModel()->contact_sequence->getDT());
    Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
    Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());
    solver_interface.GetSolution(new_times, new_states, new_inputs);
    

    // solution_t new_sol = solution_t(new_times);
    // solver_interface.getTrajectoryOptimizer()->getSolution(new_sol);

    // Eigen::MatrixXd subMatrix = new_sol.state_result.block(solver_interface.states()->q_index, 0, solver_interface.states()->nq, new_sol.state_result.cols());
    // MeshcatInterface meshcat("../examples/visualization/solution_data/");
    // meshcat.WriteTimes(new_times, "sol_times.csv");
    // meshcat.WriteJointPositions(subMatrix, "sol_states.csv");
    // meshcat.WriteMetadata(go1_location, q0_vec, "metadata.csv");
    // std::cout << "Getting constraint violations" << std::endl;
    // auto cons = solver_interface.getTrajectoryOptimizer()->getConstraintViolations(new_sol);

    // // Collect the data specific to each end effector
    // std::vector<tuple_size_t> wrench_indices;
    // std::vector<std::vector<std::string>> wrench_legend_names;
    // std::vector<std::string> ee_plot_names;
    // for (auto ee : solver_interface.getRobotModel()->getEndEffectors())
    // {
    //     wrench_indices.push_back(solver_interface.states()->frame_id_to_index_range[ee.second->frame_id]);
    //     if (ee.second->is_6d)
    //         wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}", "\\tau_{x}", "\\tau_{y}", "\\tau_{z}"});
    //     else
    //     {
    //         wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}"});
    //     }
    //     ee_plot_names.push_back("Contact Wrench of " + ee.second->frame_name);
    // }

    // GNUPlotInterface plotter(new_sol, cons);
    // plotter.PlotSolution({std::make_tuple(solver_interface.states()->q_index, solver_interface.states()->q_index + 3), std::make_tuple(solver_interface.states()->q_index + 3, solver_interface.states()->q_index + solver_interface.states()->nqb)},
    //                      wrench_indices,
    //                      {"Positions", "Orientations"},
    //                      {{"x", "y", "z"}, {"qx", "qy", "qz", "qw"}},
    //                      ee_plot_names,
    //                      wrench_legend_names);
    // plotter.PlotConstraints();

    return 0;
}