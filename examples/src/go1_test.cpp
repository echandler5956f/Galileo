#include "go1_test.h"

int main(int argc, char **argv)
{
    ConfigVector q0_vec = (ConfigVector(19) << 0., 0., 0.339, 0., 0., 0., 1., 0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3).finished();

    std::vector<int> knot_num = {20, 20, 20, 20, 20};
    std::vector<double> knot_time = {0.075, 0.45, 0.075, 0.45, 0.075};
    std::vector<uint> mask_vec = {0b1111, 0b1001, 0b1111, 0b0110};

    // std::vector<int> knot_num = {20};
    // std::vector<double> knot_time = {1.0};
    // std::vector<uint> mask_vec = {0b1111};

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
    solver_interface.LoadParameters(" ");

    casadi::DM X0 = casadi::DM::zeros(solver_interface.states()->nx, 1);
    int q0_idx = solver_interface.states()->nh + solver_interface.states()->ndh;
    for (int j = 0; j < solver_interface.states()->nq; j++)
    {
        X0(q0_idx + j) = q0_vec(j);
    }

    solver_interface.surfaces()->push_back(environment::createInfiniteGround());

    solver_interface.setContactSequence(
        knot_num,
        knot_time,
        mask_vec,
        contact_surfaces);

    solver_interface.Initialize(X0, X0);
    casadi::MXVector solution;
    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., solver_interface.robot_->contact_sequence->getDT());

    solver_interface.Update(X0, X0, solution, new_times);

    Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
    Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());
    solver_interface.solution_interface_->GetSolution(new_times, new_states, new_inputs);

    solution::solution_t new_sol = solution::solution_t(new_times, new_states, new_inputs);

    Eigen::MatrixXd subMatrix = new_sol.state_result.block(solver_interface.states()->q_index, 0, solver_interface.states()->nq, new_sol.state_result.cols());
    MeshcatInterface meshcat("../examples/visualization/solution_data/");
    meshcat.WriteTimes(new_times, "sol_times.csv");
    meshcat.WriteJointPositions(subMatrix, "sol_states.csv");
    meshcat.WriteMetadata(go1_location, q0_vec, "metadata.csv");
    std::cout << "Getting constraint violations" << std::endl;
    std::vector<std::vector<solution::constraint_evaluations_t>> cons = solver_interface.trajectory_opt_->getConstraintViolations(new_sol);

    Eigen::MatrixXd eval;
    Eigen::MatrixXd times;
    Eigen::MatrixXd lb;
    Eigen::MatrixXd ub;
    for (auto &cons_in_phase : cons)
    {
        // For each constraint type in the phase
        for (auto &constraint_type : cons_in_phase)
        {
            for (size_t i = 0; i < constraint_type.metadata.plot_titles.size(); ++i)
            {
                std::cout << constraint_type.metadata.plot_titles[i] << std::endl;
                if (constraint_type.metadata.plot_titles[i] == "Swing Velocity Constraint Vz for RL_foot")
                {
                    size_t s1 = std::get<0>(constraint_type.metadata.plot_groupings[i]);
                    size_t s2 = std::get<1>(constraint_type.metadata.plot_groupings[i]);
                    opt::solution::constraint_evaluations_t new_constraint;

                    new_constraint.evaluation = constraint_type.evaluation.block(0, s1, constraint_type.evaluation.rows(), s2 - s1);
                    new_constraint.lower_bounds = constraint_type.lower_bounds.block(0, s1, constraint_type.lower_bounds.rows(), s2 - s1);
                    new_constraint.upper_bounds = constraint_type.upper_bounds.block(0, s1, constraint_type.upper_bounds.rows(), s2 - s1);

                    eval = new_constraint.evaluation;
                    times = constraint_type.times;
                    lb = new_constraint.lower_bounds;
                    ub = new_constraint.upper_bounds;
                }
            }
        }
    }

    // save eigen matrices to a file
    std::ofstream file_eval("eval.txt");
    file_eval << eval.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
    file_eval.close();

    std::ofstream file_times("times.txt");
    file_times << times.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
    file_times.close();

    std::ofstream file_lb("lb.txt");
    file_lb << lb.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
    file_lb.close();

    std::ofstream file_ub("ub.txt");
    file_ub << ub.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
    file_ub.close();

    // Collect the data specific to each end effector
    std::vector<tuple_size_t> wrench_indices;
    std::vector<std::vector<std::string>> wrench_legend_names;
    std::vector<std::string> ee_plot_names;
    for (auto ee : solver_interface.robot_->getEndEffectors())
    {
        wrench_indices.push_back(solver_interface.states()->frame_id_to_index_range[ee.second->frame_id]);
        if (ee.second->is_6d)
            wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}", "\\tau_{x}", "\\tau_{y}", "\\tau_{z}"});
        else
        {
            wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}"});
        }
        ee_plot_names.push_back("Contact Wrench of " + ee.second->frame_name);
    }

    GNUPlotInterface plotter(new_sol, cons);
    plotter.PlotSolution({std::make_tuple(solver_interface.states()->q_index, solver_interface.states()->q_index + 3), std::make_tuple(solver_interface.states()->q_index + 3, solver_interface.states()->q_index + solver_interface.states()->nqb),
                          std::make_tuple(solver_interface.states()->q_index + solver_interface.states()->nqb, solver_interface.states()->q_index + solver_interface.states()->nqb + 3),
                          std::make_tuple(solver_interface.states()->q_index + solver_interface.states()->nqb + 3, solver_interface.states()->q_index + solver_interface.states()->nqb + 6),
                          std::make_tuple(solver_interface.states()->q_index + solver_interface.states()->nqb + 6, solver_interface.states()->q_index + solver_interface.states()->nqb + 9),
                          std::make_tuple(solver_interface.states()->q_index + solver_interface.states()->nqb + 9, solver_interface.states()->q_index + solver_interface.states()->nqb + 12)},
                         wrench_indices,
                         {"Positions", "Orientations", "FL Joints", "FR Joints", "RL Joints", "RR Joints"},
                         {{"x", "y", "z"}, {"qx", "qy", "qz", "qw"}, {"HAA", "HFE", "KFE"}, {"HAA", "HFE", "KFE"}, {"HAA", "HFE", "KFE"}, {"HAA", "HFE", "KFE"}},
                         ee_plot_names,
                         wrench_legend_names);
    plotter.PlotConstraints();

    return 0;
}