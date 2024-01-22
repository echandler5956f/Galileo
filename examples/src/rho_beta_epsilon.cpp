#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/opt/PseudospectralSegment.h"
#include <Eigen/Dense>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include "third-party/gnuplot-iostream/gnuplot-iostream.h"
#include "galileo/opt/BasicStates.h"

using namespace galileo;

int main()
{

    std::shared_ptr<opt::BasicStates> si = std::make_shared<opt::BasicStates>(std::vector<int>{4, 2});

    // Declare model variables
    casadi::SX dt = SX::sym("dt", 1);
    casadi::SX dx = SX::sym("dx", 4);

    casadi::SX x = SX::sym("x", 4);
    casadi::SX u = SX::sym("u", 2);

    double m0_p = 1.0;
    double m1_p = 1.0;
    double c0_p = 0.2;
    double c1_p = 0.2;
    double d_p = 0.1;

    casadi::SX A = horzcat(vertcat(SXVector{0., 0., (-c0_p - c1_p) / m0_p, c1_p / m1_p}),
                           vertcat(SXVector{0., 0., c1_p / m0_p, -c1_p / m1_p}),
                           vertcat(SXVector{1., 0., -2. * d_p / m0_p, d_p / m1_p}),
                           vertcat(SXVector{0., 1., d_p / m0_p, -d_p / m1_p}));

    casadi::SX B = horzcat(vertcat(SXVector{0., 0., 1. / m0_p, 0}),
                           vertcat(SXVector{0., 0., 0., 1. / m1_p}));

    // Model equations
    casadi::SX xdot = mtimes(A, x) + mtimes(B, u);

    // Objective term
    casadi::SX l = sumsqr(x - vertcat(SXVector{1., 2., 0., 0.})) + 1e-3 * sumsqr(u);

    casadi::SX x1 = SX::sym("x1", 4);
    casadi::SX x2 = SX::sym("x2", 4);

    // Continuous time dynamics (all state variables are Euclidean so we do not need to worry about these manifold operators)
    casadi::Function Fint("Fint", {x, dx, dt}, {dx});
    casadi::Function Fdif("Fdif", {x1, x2, dt}, {x2});
    casadi::Function F("F", {x, u}, {xdot});
    casadi::Function L("L", {x, u}, {l});
    casadi::Function Phi("Phi", {x}, {1e8 * sumsqr(x - vertcat(SXVector{1., 2., 0., 0.}))});

    casadi::Dict opts;
    /*DO NOT USE EXPAND OR IT WILL UNROLL ALL THE MAPS AND FOLDS*/
    // opts["ipopt.tol"] = 1e-5;
    // opts["ipopt.max_iter"] = 100;
    // opts["ipopt.print_level"] = 5;
    opts["ipopt.linear_solver"] = "ma97";
    std::shared_ptr<opt::GeneralProblemData> problem = std::make_shared<opt::GeneralProblemData>(Fint, Fdif, F, L, Phi);
    opt::TrajectoryOpt<opt::GeneralProblemData, opt::PseudospectralSegment<opt::GeneralProblemData>> traj(opts, si, problem);

    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    X0(0, 0) = 0.1;
    X0(1, 0) = 3.0;
    int d = 3;
    int N = 20;
    double T = 5.0;
    double h = T / N;

    traj.init_finite_elements(d, X0);
    auto sol = traj.optimize();

    /*TODO: Make a function to do automate the plotting within the actual library. Allow the user to include plotting functionality as a cmake option*/
    std::vector<double> times = traj.get_global_times();
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
    std::vector<double> x3_all_sol = casadi::MX::evalf(solx(2, casadi::Slice(0, solx.columns()))).get_elements();
    std::vector<double> x4_all_sol = casadi::MX::evalf(solx(3, casadi::Slice(0, solx.columns()))).get_elements();
    std::vector<double> u1_all_sol = casadi::MX::evalf(solu(0, casadi::Slice(0, solu.columns()))).get_elements();
    std::vector<double> u2_all_sol = casadi::MX::evalf(solu(1, casadi::Slice(0, solu.columns()))).get_elements();

    // std::cout << "times: " << times.size() << std::endl;
    // std::cout << "u_times: " << u_times.size() << std::endl;
    // std::cout << "x1_all_sol: " << x1_all_sol.size() << std::endl;
    // std::cout << "x2_all_sol: " << x2_all_sol.size() << std::endl;
    // std::cout << "u1_all_sol: " << u1_all_sol.size() << std::endl;
    // std::cout << "u2_all_sol: " << u2_all_sol.size() << std::endl;

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
    // gp << "plot '-' with linespoints linestyle 1 title 'x1', '-' with linespoints linestyle 2 title 'x2', '-' with linespoints linestyle 3 title 'x1_dot', '-' with linespoints linestyle 4 title 'x2_dot', '-' with linespoints linestyle 5 title 'u1', '-' with linespoints linestyle 6 title 'u2' \n";
    gp << "plot '-' with linespoints linestyle 1 title 'x1', '-' with linespoints linestyle 2 title 'x2', '-' with linespoints linestyle 3 title 'u1', '-' with linespoints linestyle 4 title 'u2' \n";
    gp.send1d(std::make_tuple(times, x1_all_sol));
    gp.send1d(std::make_tuple(times, x2_all_sol));
    // gp.send1d(std::make_tuple(times, x3_all_sol));
    // gp.send1d(std::make_tuple(times, x4_all_sol));
    gp.send1d(std::make_tuple(u_times, u1_all_sol));
    gp.send1d(std::make_tuple(u_times, u2_all_sol));

    std::cout << "Final value of x1: " << x1_all_sol[x1_all_sol.size()-1] << std::endl;
    std::cout << "Final value of x2: " << x2_all_sol[x2_all_sol.size()-1] << std::endl;
    std::cout << "Final value of x3: " << x3_all_sol[x3_all_sol.size()-1] << std::endl;
    std::cout << "Final value of x4: " << x4_all_sol[x4_all_sol.size()-1] << std::endl;
    std::cout << "Final value of u1: " << u1_all_sol[u1_all_sol.size()-1] << std::endl;
    std::cout << "Final value of u2: " << u1_all_sol[u1_all_sol.size()-1] << std::endl;

    // Eigen::MatrixXd times_mat = Eigen::Map<Eigen::MatrixXd>(times.data(), d + 1, N);
    // Eigen::MatrixXd xall_opt1 = Eigen::Map<Eigen::MatrixXd>(x1_all_sol.data(), d + 1, N);
    // Eigen::MatrixXd xall_opt2 = Eigen::Map<Eigen::MatrixXd>(x2_all_sol.data(), d + 1, N);
    // Eigen::MatrixXd uall_opt1 = Eigen::Map<Eigen::MatrixXd>(u1_all_sol.data(), d, N);
    // Eigen::MatrixXd uall_opt2 = Eigen::Map<Eigen::MatrixXd>(u2_all_sol.data(), d, N);

    // std::vector<boost::math::barycentric_rational<double>> interpolators_x1;
    // std::vector<boost::math::barycentric_rational<double>> interpolators_x2;
    // std::vector<boost::math::barycentric_rational<double>> interpolators_u1;
    // std::vector<boost::math::barycentric_rational<double>> interpolators_u2;
    // for (int k = 0; k < N; ++k)
    // {
    //     Eigen::VectorXd t_segment = times_mat.col(k);
    //     Eigen::VectorXd x1_segment = xall_opt1.col(k);
    //     Eigen::VectorXd x2_segment = xall_opt2.col(k);
    //     Eigen::VectorXd u1_segment = uall_opt1.col(k);
    //     Eigen::VectorXd u2_segment = uall_opt2.col(k);
    //     interpolators_x1.push_back(boost::math::barycentric_rational<double>(t_segment.data(), x1_segment.data(), d + 1));
    //     interpolators_x2.push_back(boost::math::barycentric_rational<double>(t_segment.data(), x2_segment.data(), d + 1));
    //     // interpolators_u1.push_back(boost::math::barycentric_rational<double>(t_segment.data(), u1_segment.data(), d + 1));
    //     // interpolators_u2.push_back(boost::math::barycentric_rational<double>(t_segment.data(), u2_segment.data(), d + 1));
    // }

    // Eigen::VectorXd t_dense = Eigen::VectorXd::LinSpaced(2 * N * (d + 1), 0, T);
    // Eigen::VectorXd x_interp_values_x1(2 * N * (d + 1) + N);
    // Eigen::VectorXd x_interp_values_x2(2 * N * (d + 1) + N);
    // Eigen::VectorXd x_interp_values_u1(2 * N * (d + 1) + N);
    // Eigen::VectorXd x_interp_values_u2(2 * N * (d + 1) + N);

    // for (int k = 0; k < N; ++k)
    // {
    //     Eigen::VectorXd t_segment = t_dense.segment(k * (2 * (d + 1)), 2 * (d + 1));
    //     x_interp_values_x1.segment(k * (2 * (d + 1)), 2 * (d + 1)) = t_segment.unaryExpr(interpolators_x1[k]);
    //     x_interp_values_x2.segment(k * (2 * (d + 1)), 2 * (d + 1)) = t_segment.unaryExpr(interpolators_x2[k]);
    //     // x_interp_values_u1.segment(k * (2 * (d + 1)), 2 * (d + 1)) = t_segment.unaryExpr(interpolators_u1[k]);
    //     // x_interp_values_u2.segment(k * (2 * (d + 1)), 2 * (d + 1)) = t_segment.unaryExpr(interpolators_u2[k]);  
    // }

    // t_dense = Eigen::VectorXd::LinSpaced(2 * N * (d + 1) + N, 0, T);

    // std::cout << t_dense.size() << std::endl;
    // std::cout << x_interp_values_x1.size() << std::endl;
    // std::cout << x_interp_values_x2.size() << std::endl;

    // // gp.send1d(std::make_tuple(t_dense, x_interp_values_x1));
    // // gp.send1d(std::make_tuple(t_dense, x_interp_values_x2));
    // // // gp.send1d(std::make_tuple(t_dense, x_interp_values_u1));
    // // // gp.send1d(std::make_tuple(t_dense, x_interp_values_u2));
    // // gp.send1d(std::make_tuple(u_times, u1_all_sol));
    // // gp.send1d(std::make_tuple(u_times, u2_all_sol));

    return 0;
}