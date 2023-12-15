#include "galileo/model/LeggedBody.h"
#include "galileo/variables/TrajectoryOpt.h"

using namespace galileo;

int main()
{

    std::shared_ptr<variables::States> si = std::make_shared<variables::States>();

    // Declare model variables
    casadi::SX dt = SX::sym("dt", 1);
    casadi::SX dx = SX::sym("dx", 2);

    casadi::SX x = SX::sym("x", 2);
    casadi::SX u = SX::sym("u");

    // Model equations
    casadi::SX xdot = vertcat((1 - pow(x(1), 2)) * x(0) - x(1) + u, x(0));

    // Objective term
    casadi::SX l = pow(x(0), 2) + pow(x(1), 2) + pow(u, 2);

    // Continuous time dynamics
    casadi::Function Fint("Fint", {x, dx, dt}, {dx});
    casadi::Function F("F", {x, u}, {xdot});
    casadi::Function L("L", {x, u}, {l});
    casadi::Function Phi("Phi", {x}, {0});

    casadi::Dict opts;
    opts["ipopt.tol"] = 1e-5;
    // opts["ipopt.max_iter"] = 100;
    // opts["ipopt.print_level"] = 5;
    opts["expand"] = true;
    opts["ipopt.linear_solver"] = "ma97";
    std::shared_ptr<variables::ProblemData> problem = std::make_shared<variables::ProblemData>(Fint, F, L, Phi);
    variables::TrajectoryOpt traj(opts, si, problem);

    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    X0(1,0) = 1;

    traj.init_finite_elements(9, X0);

    auto sol = traj.optimize();
    std::cout << sol["x"] << std::endl;

    return 0;
}