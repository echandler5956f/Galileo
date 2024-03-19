#include "simple_test.h"

int main()
{
    std::shared_ptr<BasicStates> si = std::make_shared<BasicStates>(2, 1);

    // Declare model variables
    casadi::SX dt = casadi::SX::sym("dt", 1);
    casadi::SX dx = casadi::SX::sym("dx", 2, 1);

    casadi::SX x = casadi::SX::sym("x", 2, 1);
    casadi::SX u = casadi::SX::sym("u", 1);

    // Model equations
    casadi::SX xdot = vertcat((1 - pow(x(1), 2)) * x(0) - x(1) + u, x(0));

    // Objective term
    casadi::SX l = pow(x(0), 2) + pow(x(1), 2) + pow(u, 2);

    casadi::SX x1 = casadi::SX::sym("x1", 2, 1);
    casadi::SX x2 = casadi::SX::sym("x2", 2, 1);

    // Continuous time dynamics (all state variables are Euclidean so we do not need to worry about these manifold operators)
    casadi::Function Fint("Fint", {x, dx, dt}, {dx});
    casadi::Function Fdif("Fdif", {x1, x2, dt}, {x2});
    casadi::Function F("F", {x, u}, {xdot});
    casadi::Function L("L", {x, u}, {l});
    casadi::Function Phi("Phi", {x}, {0});

    std::shared_ptr<GeneralProblemData> gp_data = std::make_shared<GeneralProblemData>(Fint, Fdif, L, Phi);
    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    X0(1, 0) = 1;

    int d = 9;
    int N = 20;
    double T = 10.;
    double h = T / N;

    std::shared_ptr<BasicSequence> seq = std::make_shared<BasicSequence>();
    BasicMode mode;
    seq->addPhase(mode, N, T, F);

    std::shared_ptr<RosenbrockProblemData> problem = std::make_shared<RosenbrockProblemData>(gp_data, seq, si, x, u, dt);
    // std::shared_ptr<ConstraintBuilder<RosenbrockProblemData>> simple_builder =
    //     std::make_shared<SimpleConstraintBuilder<RosenbrockProblemData>>();

    std::vector<std::shared_ptr<ConstraintBuilder<RosenbrockProblemData>>> builders = {};
    std::shared_ptr<DecisionDataBuilder<RosenbrockProblemData>> decision_builder = std::make_shared<SimpleDecisionDataBuilder<RosenbrockProblemData>>();

    casadi::Dict opts;
    opts["ipopt.linear_solver"] = "ma97";
    opts["ipopt.ma97_order"] = "metis";
    opts["ipopt.fixed_variable_treatment"] = "make_constraint";
    opts["ipopt.max_iter"] = 250;
    TrajectoryOpt<RosenbrockProblemData, BasicMode> traj(problem, seq, builders, decision_builder, opts, "ipopt");

    traj.initFiniteElements(d, X0);
    casadi::MXVector sol = traj.optimize();

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(20, 0., seq->getDT());
    std::vector<opt::solution::solution_segment_data_t> solution_segments_;
    std::shared_ptr<opt::solution::Solution> solution_interface_ = std::make_shared<opt::solution::Solution>();
    solution_segments_.clear();
    traj.getSolutionSegments(solution_segments_, new_times);
    solution_interface_->UpdateSolution(solution_segments_);
    Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(si->nx, new_times.size());
    Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(si->nu, new_times.size());
    solution_interface_->GetSolution(new_times, new_states, new_inputs);
    solution::solution_t new_sol = solution::solution_t(new_times, new_states, new_inputs);
    auto cons = traj.getConstraintViolations(new_sol);

    // GNUPlotInterface plotter(new_sol, cons);
    // plotter.PlotSolution({std::make_tuple(0, si->nx)},
    //                      {std::make_tuple(0, si->nu)},
    //                      {"States"},
    //                      {{"x1", "x2"}},
    //                      {"Input"},
    //                      {{"u"}});
    // plotter.PlotConstraints();

    Gnuplot gp;
    gp << "set terminal pngcairo\n";
    gp << "set output 'Van der Pol.png'\n";
    gp << "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#FFFFFF' behind\n";
    gp << "set xlabel 'Time'\n";
    gp << "set xrange [" << new_sol.times(0) << ":" << new_sol.times(new_sol.times.size() - 1) << "]\n";
    gp << "set ylabel 'Values'\n";
    gp << "set title 'Van der Pol Oscillator'\n";
    gp << "plot '-' with lines linetype 1 linewidth 2 linecolor rgb '#ACC18A' title 'x_1', "
       << "'-' with lines linetype 1 linewidth 2 linecolor rgb '#7D8CC4' title 'x_2', "
       << "'-' with lines linetype 1 dashtype 2 linecolor rgb '#B8336A' title 'u'\n";
    Eigen::MatrixXd new_sol_state_transpose = new_sol.state_result.transpose();
    std::vector<double> std_times_vector(new_sol.times.data(), new_sol.times.data() + new_sol.times.size());

    for (auto j = 0; j < new_sol_state_transpose.cols(); j++)
    {
        Eigen::VectorXd colmatrix_state = new_sol_state_transpose.col(j).matrix();

        std::vector<double> std_state_vector(colmatrix_state.data(), colmatrix_state.data() + colmatrix_state.size());
        gp.send1d(std::make_tuple(std_times_vector, std_state_vector));
    }

    std::vector<double> std_input_vector(new_sol.input_result.data(), new_sol.input_result.data() + new_sol.input_result.size());
    gp.send1d(std::make_tuple(std_times_vector, std_input_vector));

    return 0;
}