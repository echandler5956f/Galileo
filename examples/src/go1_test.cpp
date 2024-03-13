#include "go1_test.h"

int main(int argc, char **argv)
{

    ConfigVector q0_vec = (ConfigVector(19) << 0., 0., 0.339, 0., 0., 0., 1., 0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3).finished();

    std::vector<int> knot_num = {20, 20, 20, 20, 20};                      //, 20, 20, 20, 20};
    std::vector<double> knot_time = {0.05, 0.3, 0.05, 0.3, 0.05};          // , 0.3, 0.05, 0.3, 0.05};
    std::vector<uint> mask_vec = {0b1111, 0b1101, 0b1111, 0b1011, 0b1111}; //, 0b1110, 0b1111, 0b0111, 0b1111}; // static walk

    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;

    for (int i = 0; i < mask_vec.size(); i++)
    {
        std::vector<galileo::legged::environment::SurfaceID> tmp;
        for (int j = 0; j < 4; j++)
        {
            if (mask_vec[i] & (1 << j))
            {
                tmp.push_back(0);
            }
            else
            {
                tmp.push_back(galileo::legged::environment::NO_SURFACE);
            }
        }
        contact_surfaces.push_back(tmp);
    }

    galileo::legged::LeggedInterface solver_interface;

    solver_interface.LoadModel(go1_location, end_effector_names);
    solver_interface.LoadParameters(" ");

    casadi::DM X0 = casadi::DM::zeros(solver_interface.states()->nx, 1);
    int j = solver_interface.states()->nh + solver_interface.states()->ndh;
    int q0_idx = solver_interface.states()->nh + solver_interface.states()->ndh;
    for (int j = 0; j < solver_interface.states()->nq; j++)
    {
        X0(q0_idx + j) = q0_vec(j);
        ++j;
    }

    solver_interface.surfaces()->push_back(environment::createInfiniteGround());

    solver_interface.setContactSequence(
        knot_num,
        knot_time,
        mask_vec,
        contact_surfaces);

    solver_interface.Initialize(X0, X0);
    casadi::MXVector solution;
    solver_interface.Update(X0, X0, solution);

    // std::cout << "Total duration: " << robot.contact_sequence->getDT() << std::endl;
    // Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., robot.contact_sequence->getDT());
    // std::vector<double> tmp = traj.getGlobalTimes();
    // double threshold = 1e-6; // Adjust this value as needed

    // // Custom comparison function
    // auto closeEnough = [threshold](double a, double b) {
    //     return std::abs(a - b) < threshold;
    // };

    // std::sort(tmp.begin(), tmp.end());
    // auto last = std::unique(tmp.begin(), tmp.end(), closeEnough);
    // tmp.erase(last, tmp.end());

    // Eigen::VectorXd new_times = Eigen::Map<Eigen::VectorXd>(tmp.data(), tmp.size());
    // std::cout << "New times: " << new_times << std::endl;

    // solution_t new_sol = solution_t(new_times);
    // traj.getSolution(new_sol);

    // Eigen::MatrixXd subMatrix = new_sol.state_result.block(si->nh + si->ndh, 0, si->nq, new_sol.state_result.cols());
    // MeshcatInterface meshcat("../examples/visualization/");
    // meshcat.WriteTimes(new_times, "sol_times.csv");
    // meshcat.WriteJointPositions(subMatrix, "sol_states.csv");
    // meshcat.WriteMetadata(go1_location, q0_vec, "metadata.csv");
    // std::cout << "Getting constraint violations" << std::endl;
    // auto cons = traj.getConstraintViolations(new_sol);

    // // Collect the data specific to each end effector
    // std::vector<tuple_size_t> wrench_indices;
    // std::vector<std::vector<std::string>> wrench_legend_names;
    // std::vector<std::string> ee_plot_names;
    // for (auto ee : robot.getEndEffectors())
    // {
    //     // auto range = si->frame_id_to_index_range[ee.second->frame_id];
    //     // wrench_indices.push_back(std::make_tuple(std::get<1>(range) - 1, std::get<1>(range)));
    //     wrench_indices.push_back(si->frame_id_to_index_range[ee.second->frame_id]);
    //     if (ee.second->is_6d)
    //         wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}", "\\tau_{x}", "\\tau_{y}", "\\tau_{z}"});
    //     else
    //     {
    //         wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}"});
    //         // wrench_legend_names.push_back({"F_{z}"});
    //     }
    //     ee_plot_names.push_back("Contact Wrench of " + ee.second->frame_name);
    // }

    // GNUPlotInterface plotter(new_sol, cons);
    // plotter.PlotSolution({std::make_tuple(si->nh + si->ndh, si->nh + si->ndh + 3), std::make_tuple(si->nh + si->ndh + 3, si->nh + si->ndh + si->nqb)},
    //                      wrench_indices,
    //                      {"Positions", "Orientations"},
    //                      {{"x", "y", "z"}, {"qx", "qy", "qz", "qw"}},
    //                      ee_plot_names,
    //                      wrench_legend_names);
    // plotter.PlotConstraints();

    // Eigen::MatrixXd new_state(new_sol.state_result.rows() - 1, new_sol.state_result.cols());
    // Eigen::MatrixXd euler_angles = galileo::math::quaternion2Euler(new_sol.state_result.block(si->nh + si->ndh + 3, 0, 4, new_sol.state_result.cols()), galileo::math::zyx);
    // new_state << new_sol.state_result.topRows(si->nh + si->ndh + 3), euler_angles, new_sol.state_result.bottomRows(si->nx - (si->nh + si->ndh + si->nqb));
    // std::cout << "new_state: " << new_state << std::endl;

    return 0;
}