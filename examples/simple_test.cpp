#include "galileo/model/LeggedBody.h"
#include "galileo/opt/TrajectoryOpt.h"
#include <Eigen/Dense>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include "third-party/gnuplot-iostream/gnuplot-iostream.h"

using namespace galileo;

int main()
{

    std::shared_ptr<opt::States> si = std::make_shared<opt::States>();

    // Declare model variables
    casadi::SX dt = SX::sym("dt", 1);
    casadi::SX dx = SX::sym("dx", 2);

    casadi::SX x = SX::sym("x", 2);
    casadi::SX u = SX::sym("u");

    // Model equations
    casadi::SX xdot = vertcat((1 - pow(x(1), 2)) * x(0) - x(1) + u, x(0));

    // Objective term
    casadi::SX l = pow(x(0), 2) + pow(x(1), 2) + pow(u, 2);

    casadi::SX x1 = SX::sym("x1", 2);
    casadi::SX x2 = SX::sym("x2", 2);

    // Continuous time dynamics (all state variables are Euclidean so we do not need to worry about these manifold operators)
    casadi::Function Fint("Fint", {x, dx, dt}, {dx});
    casadi::Function Fdif("Fdif", {x1, x2, dt}, {x2});
    casadi::Function F("F", {x, u}, {xdot});
    casadi::Function L("L", {x, u}, {l});
    casadi::Function Phi("Phi", {x}, {0});

    casadi::Dict opts;
    /*DO NOT USE EXPAND OR IT WILL UNROLL ALL THE MAPS AND FOLDS*/
    // opts["ipopt.tol"] = 1e-5;
    // opts["ipopt.max_iter"] = 100;
    // opts["ipopt.print_level"] = 5;
    opts["ipopt.linear_solver"] = "ma97";
    std::shared_ptr<opt::GeneralProblemData> problem = std::make_shared<opt::GeneralProblemData>(Fint, Fdif, F, L, Phi);
    opt::TrajectoryOpt traj(opts, si, problem);

    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    X0(1, 0) = 1;
    int d = 3;
    int N = 20;
    double T = 10.0;
    double h = T / N;

    traj.init_finite_elements(d, X0);
    auto sol = traj.optimize();

    /*TODO: Make a function to do automate the plotting within the actual library. Allow the user to include plotting functionality as a cmake option*/
    std::vector<double> times = traj.get_all_times();
    std::vector<double> u_times;
    for (int i = 0; i < times.size(); ++i)
    {
        if (i % (d + 1) != d) // Skip every dth element
        {
            u_times.push_back(times[i]);
        }
    }
    // std::cout << "sol: " << sol << std::endl;
    auto solx = sol[0];
    auto solu = sol[1];

    std::vector<double> x1_all_sol = casadi::MX::evalf(solx(0, casadi::Slice(0, solx.columns()))).get_elements();
    std::vector<double> x2_all_sol = casadi::MX::evalf(solx(1, casadi::Slice(0, solx.columns()))).get_elements();
    std::vector<double> u_all_sol = casadi::MX::evalf(solu(0, casadi::Slice(0, solu.columns()))).get_elements();

    std::cout << "times: " << times.size() << std::endl;
    std::cout << "u_times: " << u_times.size() << std::endl;
    std::cout << "x1_all_sol: " << x1_all_sol.size() << std::endl;
    std::cout << "x2_all_sol: " << x2_all_sol.size() << std::endl;
    std::cout << "u_all_sol: " << u_all_sol.size() << std::endl;

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
    gp.send1d(std::make_tuple(times, x1_all_sol));
    gp.send1d(std::make_tuple(times, x2_all_sol));
    gp.send1d(std::make_tuple(u_times, u_all_sol));

    // // Something is buggy with this interpolation. Do not trust it too much
    // Eigen::MatrixXd xall_opt1 = Eigen::Map<Eigen::MatrixXd>(x1_all_sol.data(), d + 1, N);
    // Eigen::MatrixXd xall_opt2 = Eigen::Map<Eigen::MatrixXd>(x2_all_sol.data(), d + 1, N);
    // Eigen::MatrixXd uall_opt = Eigen::Map<Eigen::MatrixXd>(u_all_sol.data(), d, N);

    // std::vector<boost::math::barycentric_rational<double>> interpolators_x1;
    // std::vector<boost::math::barycentric_rational<double>> interpolators_x2;
    // std::vector<boost::math::barycentric_rational<double>> interpolators_u;
    // for (int k = 0; k < N; ++k)
    // {
    //     Eigen::VectorXd t_segment = Eigen::VectorXd::LinSpaced(d + 1, k * h, (k + 1) * h);
    //     Eigen::VectorXd x1_segment = xall_opt1.col(k);
    //     Eigen::VectorXd x2_segment = xall_opt2.col(k);
    //     Eigen::VectorXd u_segment = uall_opt.col(k);

    //     interpolators_x1.push_back(boost::math::barycentric_rational<double>(t_segment.data(), x1_segment.data(), d + 1));
    //     interpolators_x2.push_back(boost::math::barycentric_rational<double>(t_segment.data(), x2_segment.data(), d + 1));
    //     interpolators_u.push_back(boost::math::barycentric_rational<double>(t_segment.data(), u_segment.data(), d));
    // }

    // Eigen::VectorXd t_dense = Eigen::VectorXd::LinSpaced(2 * N * (d + 1), 0, T);
    // Eigen::VectorXd x_interp_values_x1(2 * N * (d + 1) + N);
    // Eigen::VectorXd x_interp_values_x2(2 * N * (d + 1) + N);
    // Eigen::VectorXd x_interp_values_u(2 * N * (d + 1) + N);

    // for (int k = 0; k < N; ++k)
    // {
    //     Eigen::VectorXd t_segment = t_dense.segment(k * (2 * (d + 1)), 2 * (d + 1));
    //     x_interp_values_x1.segment(k * (2 * (d + 1)), 2 * (d + 1)) = t_segment.unaryExpr(interpolators_x1[k]);
    //     x_interp_values_x2.segment(k * (2 * (d + 1)), 2 * (d + 1)) = t_segment.unaryExpr(interpolators_x2[k]);
    //     x_interp_values_u.segment(k * (2 * (d + 1)), 2 * (d + 1)) = t_segment.unaryExpr(interpolators_u[k]);
    // }

    // t_dense = Eigen::VectorXd::LinSpaced(2 * N * (d + 1) + N, 0, T);

    // std::cout << t_dense.size() << std::endl;
    // std::cout << x_interp_values_x1.size() << std::endl;
    // std::cout << x_interp_values_x2.size() << std::endl;

    // gp.send1d(std::make_tuple(t_dense, x_interp_values_x1));
    // gp.send1d(std::make_tuple(t_dense, x_interp_values_x2));
    // gp.send1d(std::make_tuple(t_dense, x_interp_values_u));
    return 0;
}