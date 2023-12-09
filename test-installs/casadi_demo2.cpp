#include <pinocchio/autodiff/casadi.hpp>
#include <iostream>
#include <vector>
#include <string>

using namespace casadi;
using namespace std;

int main(int argc, char **argv)
{
  // Degree of interpolating polynomial
  int d = 2;

  // if (argc > 1) {
  //   d = std::atoi(argv[1]);
  //   std::cout << "Used arg" << std::endl;
  // }
  // else {
  //   std::cout << "Did not use arg" << std::endl;
  // }

  // Choose collocation points
  auto troot = collocation_points(d, "radau");
  troot.insert(troot.begin(), 0);
  Eigen::VectorXd tau_root = Eigen::Map<Eigen::VectorXd>(troot.data(), troot.size());

  // Coefficients of the collocation equation
  Eigen::MatrixXd C(d + 1, d + 1);

  // Coefficients of the continuity equation
  Eigen::VectorXd D(d + 1);

  // Coefficients of the quadrature function
  Eigen::VectorXd B(d + 1);

  // For all collocation points
  for (int j = 0; j < d + 1; ++j)
  {

    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    Polynomial p = 1;
    for (int r = 0; r < d + 1; ++r)
    {
      if (r != j)
      {
        p *= Polynomial(-tau_root(r), 1) / (tau_root(j) - tau_root(r));
      }
    }
    // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D(j) = p(1.0);

    // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    Polynomial dp = p.derivative();
    for (int r = 0; r < d + 1; ++r)
    {
      C(j, r) = dp(tau_root(r));
    }
    auto pint = p.anti_derivative();
    B(j) = pint(1.0);
  }

  // Time horizon
  double T = 10.0;

  // Declare model variables
  SX x = SX::sym("x", 2);
  SX u = SX::sym("u");

  // Model equations
  SX xdot = vertcat((1 - pow(x(1), 2)) * x(0) - x(1) + u, x(0));

  // Objective term
  SX L = pow(x(0), 2) + pow(x(1), 2) + pow(u, 2);

  // Continuous time dynamics
  Function f("f", {x, u}, {xdot, L});

  // Control discretization
  int N = 2; // number of control intervals
  double h = T / N;

  // Variables for a single collocation interval
  vector<SX> Xc;
  for (int j = 0; j < d; ++j)
  {
    Xc.push_back(SX::sym("Xc_" + to_string(j), 2));
  }
  SX X0 = SX::sym("X0", 2);
  SX P = SX::sym("P");
  SX Lc = SX::sym("L");

  // State at the end of the collocation interval
  SX Xf = D(0) * X0;
  for (int j = 0; j < d; ++j)
  {
    Xf += D(j + 1) * Xc[j];
  }

  // Nonlinear equations and cost contribution for collocation interval
  vector<SX> eq;
  SX Qf = 0.0;
  for (int j = 0; j < d; ++j)
  {
    // Expression for the state derivative at the collocation point
    SX xp = C(0, j + 1) * X0;
    for (int r = 0; r < d; ++r)
    {
      xp += C(r + 1, j + 1) * Xc[r];
    }

    // Evaluate the function at the collocation points
    SX fj = f(vector<SX>{Xc[j], P}).at(0);

    // Append collocation equations
    eq.push_back(h * f(vector<SX>{Xc[j], P}).at(0) - xp);

    // Add cost contribution
    Qf += B(j + 1) * f(vector<SX>{Xc[j], P}).at(1) * h;
  }

  // Implicit discrete-time dynamics
  Function F("F", {vertcat(Xc), X0, P}, {vertcat(eq), Xf, Qf});
  Function Feq("feq", {vertcat(Xc), X0, P}, {vertcat(eq)});
  Function Fxf("feq", {vertcat(Xc), X0, P}, {Xf});
  Function Fxq("fxq", {Lc, vertcat(Xc), X0, P}, {Lc + Qf});

  // Optimization
  Opti opti;

  MX J = 0.0;
  vector<MX> var_xs(N + 1);
  vector<MX> var_us(N);
  vector<MX> var_xcs(N);

  var_xs[0] = opti.variable(2);
  opti.subject_to(var_xs[0] == vector<double>{0, 1});
  opti.set_initial(var_xs[0], vector<double>{0, 1});

  Function feqmap = Feq.map(N, "openmp");
  Function fxfmap = Fxf.map(N, "openmp");
  Function fxqfold = Fxq.fold(N);

  vector<double> lb;
  vector<double> ub;
  for (int r = 0; r < d; ++r)
  {
    lb.push_back(-numeric_limits<double>::infinity());
    lb.push_back(-0.25);
    ub.push_back(numeric_limits<double>::infinity());
    ub.push_back(numeric_limits<double>::infinity());
  }

  for (int k = 0; k < N; ++k)
  {
    var_us[k] = opti.variable(1);
    opti.bounded(-1, var_us[k], 1);

    var_xcs[k] = opti.variable(2 * d);
    opti.bounded(lb, var_xcs[k], ub);

    var_xs[k + 1] = opti.variable(2);
    opti.bounded(vector<double>{-0.25, -numeric_limits<double>::infinity()}, var_xs[k + 1], vector<double>{numeric_limits<double>::infinity(), numeric_limits<double>::infinity()});
  }

  MX xsoffset = var_xs[1];
  MX xs = var_xs[0];
  MX us = var_us[0];
  MX xcs = var_xcs[0];

  for (int k = 1; k < N; ++k)
  {
    xs = horzcat(xs, var_xs[k]);
    us = horzcat(us, var_us[k]);
    xcs = horzcat(xcs, var_xcs[k]);
    xsoffset = horzcat(xsoffset, var_xs[k + 1]);
  }
  std::cout << "feqmap: " << feqmap(vector<MX>{xcs, xs, us}).at(0).size() << std::endl;
  opti.subject_to(feqmap(vector<MX>{xcs, xs, us}).at(0) == 0);
  opti.subject_to(fxfmap(vector<MX>{xcs, xs, us}).at(0) - xsoffset == 0);
  J = fxqfold(vector<MX>{J, xcs, xs, us}).at(0);

  opti.minimize(J);
  Dict p_opts = {{"expand", true}};
  Dict s_opts = {{"max_iter", 200}, 
  // {"linear_solver", "ma97"}
  };
  opti.solver("ipopt", p_opts, s_opts);
  auto sol = opti.solve();

  return 0;
}
