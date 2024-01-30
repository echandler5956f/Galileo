#include "simple_test.h"

int main()
{

    std::shared_ptr<opt::BasicStates> si = std::make_shared<opt::BasicStates>(std::vector<int>{2, 1});

    // Declare model variables
    casadi::SX dt = casadi::SX::sym("dt", 1);
    casadi::SX dx = casadi::SX::sym("dx", 2);

    casadi::SX x = casadi::SX::sym("x", 2);
    casadi::SX u = casadi::SX::sym("u");

    // Model equations
    casadi::SX xdot = vertcat((1 - pow(x(1), 2)) * x(0) - x(1) + u, x(0));

    // Objective term
    casadi::SX l = pow(x(0), 2) + pow(x(1), 2) + pow(u, 2);

    casadi::SX x1 = casadi::SX::sym("x1", 2);
    casadi::SX x2 = casadi::SX::sym("x2", 2);

    // Continuous time dynamics (all state variables are Euclidean so we do not need to worry about these manifold operators)
    casadi::Function Fint("Fint", {x, dx, dt}, {dx});
    casadi::Function Fdif("Fdif", {x1, x2, dt}, {x2});
    casadi::Function F("F", {x, u}, {xdot});
    casadi::Function L("L", {x, u}, {l});
    casadi::Function Phi("Phi", {x}, {0});

    std::shared_ptr<opt::GeneralProblemData> gp_data = std::make_shared<opt::GeneralProblemData>(Fint, Fdif, F, L, Phi);
    std::shared_ptr<SimpleProblemData> problem = std::make_shared<SimpleProblemData>();
    problem->gp_data = gp_data;
    problem->x = x;
    problem->u = u;
    problem->states = si;
    problem->num_knots = 20;
    std::shared_ptr<opt::ConstraintData> constraint_data = std::make_shared<opt::ConstraintData>();
    std::shared_ptr<RosenbrockProblemData> rosenbrock_problem_data = std::make_shared<RosenbrockProblemData>();

    std::shared_ptr<opt::ConstraintBuilder<SimpleProblemData>> rosenbrock_builder =
        std::make_shared<RosenbrockConstraintBuilder<SimpleProblemData>>();

    std::vector<std::shared_ptr<opt::ConstraintBuilder<SimpleProblemData>>> builders = {rosenbrock_builder};

    casadi::Dict opts;
    /*DO NOT USE EXPAND OR IT WILL UNROLL ALL THE MAPS AND FOLDS*/
    // opts["ipopt.tol"] = 1e-5;
    // opts["ipopt.max_iter"] = 100;
    // opts["ipopt.print_level"] = 5;
    opts["ipopt.linear_solver"] = "ma97";
    opt::TrajectoryOpt<SimpleProblemData> traj(problem, builders, opts);

    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    X0(1, 0) = 1;
    int d = 3;
    int N = 20;
    double T = 10.0;
    double h = T / N;

    traj.init_finite_elements(d, X0);
    auto sol = traj.optimize();

    /*TODO: Make a function to do automate the plotting within the actual library. Allow the user to include plotting functionality as a cmake option*/
    std::vector<double> times = traj.get_global_times();
    std::vector<double> u_times;
    for (int i = 0; i < times.size(); ++i)
    {
        if (i % (d + 1) != d) // Skip every dth element, because u does not include knot points
        {
            u_times.push_back(times[i]);
        }
    }
    // std::cout << "sol: " << sol << std::endl;
    auto solx = sol[0];
    auto solu = sol[1];

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(100, 0, 10.0);
    // Eigen::VectorXd new_times = Eigen::Map<Eigen::VectorXd>(times.data(), times.size());
    Eigen::MatrixXd new_sol = traj.get_solution(new_times);

    std::vector<double> x1_all_sol = casadi::MX::evalf(solx(0, casadi::Slice(0, solx.columns()))).get_elements();
    std::vector<double> x2_all_sol = casadi::MX::evalf(solx(1, casadi::Slice(0, solx.columns()))).get_elements();
    std::vector<double> u_all_sol = casadi::MX::evalf(solu(0, casadi::Slice(0, solu.columns()))).get_elements();

    // std::cout << "times: " << times.size() << std::endl;
    // std::cout << "u_times: " << u_times.size() << std::endl;
    // std::cout << "x1_all_sol: " << x1_all_sol.size() << std::endl;
    // std::cout << "x2_all_sol: " << x2_all_sol.size() << std::endl;
    // std::cout << "u_all_sol: " << u_all_sol.size() << std::endl;

    // Create a Gnuplot object
    Gnuplot gp;

    // Set the style of the plot
    // gp << "set style line 1 lc rgb '#2ca02c' lt 1 lw 2 pt 7 ps 0.5\n"; // Cooked asparagus green
    // gp << "set style line 2 lc rgb '#9467bd' lt 1 lw 2 pt 7 ps 0.5\n"; // Muted purple
    // gp << "set style line 3 lc rgb '#8c564b' lt 1 lw 2 pt 7 ps 0.5\n"; // Chestnut brown

    // Set labels and title
    gp << "set xlabel 'Time'\n";
    gp << "set ylabel 'Values'\n";
    gp << "set title 'Solution Plot'\n";

    // Plot the data
    gp << "plot '-' with linespoints linestyle 1 title 'x1', '-' with linespoints linestyle 2 title 'x2', '-' with linespoints linestyle 3 title 'u'\n";
    gp.send1d(std::make_tuple(new_times, new_sol.row(0)));
    gp.send1d(std::make_tuple(new_times, new_sol.row(1)));
    // gp.send1d(std::make_tuple(times, x1_all_sol));
    // gp.send1d(std::make_tuple(times, x2_all_sol));
    // gp.send1d(std::make_tuple(u_times, u_all_sol));

    return 0;
}