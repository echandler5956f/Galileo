#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/opt/PseudospectralSegment.h"
#include <Eigen/Dense>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include "third-party/gnuplot-iostream/gnuplot-iostream.h"
#include "galileo/opt/BasicStates.h"
#include <cstdlib> // for rand() and srand()

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

    std::vector<double> total_times;
    std::vector<double> total_x1;
    std::vector<double> total_x2;
    std::vector<double> total_x3;
    std::vector<double> total_x4;
    std::vector<double> total_u1;
    std::vector<double> total_u2;

    for (int i = 0; i < 100; ++i)
    {
        casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
        X0(0, 0) = static_cast<double>(rand()) / RAND_MAX * 3.;
        X0(1, 0) = static_cast<double>(rand()) / RAND_MAX * 3.;
        X0(2, 0) = static_cast<double>(rand()) / RAND_MAX * 3.;
        X0(3, 0) = static_cast<double>(rand()) / RAND_MAX * 3.;
        int d = 9;
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

        total_times.insert(total_times.end(), times.begin(), times.end());
        total_x1.insert(total_x1.end(), x1_all_sol.begin(), x1_all_sol.end());
        total_x2.insert(total_x2.end(), x2_all_sol.begin(), x2_all_sol.end());
        total_x3.insert(total_x3.end(), x3_all_sol.begin(), x3_all_sol.end());
        total_x4.insert(total_x4.end(), x4_all_sol.begin(), x4_all_sol.end());
        total_u1.insert(total_u1.end(), u1_all_sol.begin(), u1_all_sol.end());
        total_u2.insert(total_u2.end(), u2_all_sol.begin(), u2_all_sol.end());
    }

    Eigen::MatrixXd plotsol3d_x1(2, total_times.size());
    Eigen::MatrixXd plotsol3d_x2(2, total_times.size());
    Eigen::MatrixXd statePlot1(2, total_times.size());
    Eigen::MatrixXd statePlot2(2, total_times.size());

    // Fill the matrix with data from the vectors
    for (int i = 0; i < total_times.size(); ++i)
    {
        plotsol3d_x1(0, i) = total_times[i];
        plotsol3d_x1(1, i) = total_x1[i];

        plotsol3d_x2(0, i) = total_times[i];
        plotsol3d_x2(1, i) = total_x2[i];

        statePlot1(0, i) = total_x1[i];
        statePlot1(1, i) = total_x3[i];

        statePlot2(0, i) = total_x2[i];
        statePlot2(1, i) = total_x4[i];        
    }

    // Create a Gnuplot object
    Gnuplot gp;

    // Set labels and title
    gp << "set xlabel 'time'\n";
    gp << "set ylabel 'x1'\n";
    gp << "set title 'Solution Plot'\n";

    // Plot the data
    gp << "plot '-' with points linestyle 1 title 'x1' \n";
    gp.send1d(plotsol3d_x1.transpose());

    Gnuplot gp2;

    // Set labels and title
    gp2 << "set xlabel 'time'\n";
    gp2 << "set ylabel 'x2'\n";
    gp2 << "set title 'Solution Plot'\n";

    // Plot the data
    gp2 << "plot '-' with points linestyle 2 title 'x2' \n";
    gp2.send1d(plotsol3d_x2.transpose());

    
    Gnuplot gp3;

    // Set labels and title
    gp3 << "set xlabel 'x1'\n";
    gp3 << "set ylabel 'x1_{dot}'\n";
    gp3 << "set title 'Solution Plot'\n";

    // Plot the data
    gp3 << "plot '-' with points linestyle 2 title 'x1 vs x1_{dot}' \n";
    gp3.send1d(statePlot1.transpose());

    
    Gnuplot gp4;

    // Set labels and title
    gp4 << "set xlabel 'x2'\n";
    gp4 << "set ylabel 'x2_{dot}'\n";
    gp4 << "set title 'Solution Plot'\n";

    // Plot the data
    gp4 << "plot '-' with points linestyle 2 title 'x2 vs x2_{dot}' \n";
    gp4.send1d(statePlot2.transpose());

    return 0;
}