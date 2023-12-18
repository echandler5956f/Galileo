// #include <casadi/casadi.hpp>
// #include <iostream>
// #include <vector>

// using namespace casadi;
// using namespace std;

// int main(int argc, char *argv[])
// {
//   // Degree of interpolating polynomial
//   int d = 3;

//   // Choose collocation points
//   vector<double> tau_root = collocation_points(d, "legendre");
//   tau_root.insert(tau_root.begin(), 0);

//   // Coefficients of the collocation equation
//   vector<vector<double>> C(d + 1, vector<double>(d + 1, 0));

//   // Coefficients of the continuity equation
//   vector<double> D(d + 1, 0);

//   // Coefficients of the quadrature function
//   vector<double> B(d + 1, 0);

//   // For all collocation points
//   for (int j = 0; j < d + 1; ++j)
//   {

//     // Construct Lagrange polynomials to get the polynomial basis at the collocation point
//     Polynomial p = 1;
//     for (int r = 0; r < d + 1; ++r)
//     {
//       if (r != j)
//       {
//         p *= Polynomial(-tau_root[r], 1) / (tau_root[j] - tau_root[r]);
//       }
//     }

//     // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
//     D[j] = p(1.0);

//     // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
//     Polynomial dp = p.derivative();
//     for (int r = 0; r < d + 1; ++r)
//     {
//       C[j][r] = dp(tau_root[r]);
//     }
//     auto pint = p.anti_derivative();
//     B[j] = pint(1.0);
//   }

//   // Time horizon
//   double T = 10.0;

//   // Declare model variables
//   SX x = SX::sym("x", 2);
//   SX u = SX::sym("u");

//   // Model equations
//   SX xdot = vertcat((1 - pow(x(1), 2)) * x(0) - x(1) + u, x(0));

//   // Objective term
//   SX L = pow(x(0), 2) + pow(x(1), 2) + pow(u, 2);

//   // Continuous time dynamics
//   Function f("f", {x, u}, {xdot, L});

//   // Control discretization
//   int N = 20; // number of control intervals
//   double h = T / N;

//   // Variables for a single collocation interval
//   vector<SX> Xc;
//   for (int j = 0; j < d; ++j)
//   {
//     Xc.push_back(SX::sym("Xc_" + to_string(j), 2));
//   }
//   SX X0 = SX::sym("X0", 2);
//   SX U = SX::sym("U");
//   SX Lc = SX::sym("L");

//   // State at the end of the collocation interval
//   SX Xf = D[0] * X0;
//   for (int j = 0; j < d; ++j)
//   {
//     Xf += D[j + 1] * Xc[j];
//   }

//   // Nonlinear equations and cost contribution for collocation interval
//   vector<SX> eq;
//   SX Qf = 0.0;
//   for (int j = 0; j < d; ++j)
//   {
//     // Expression for the state derivative at the collocation point
//     SX xp = C[0][j + 1] * X0;
//     for (int r = 0; r < d; ++r)
//     {
//       xp += C[r + 1][j + 1] * Xc[r];
//     }

//     // Evaluate the function at the collocation points
//     SX fj = f(vector<SX>{Xc[j], U}).at(0);

//     // Append collocation equations
//     eq.push_back(h * f(vector<SX>{Xc[j], U}).at(0) - xp);

//     // Add cost contribution
//     Qf += B[j + 1] * f(vector<SX>{Xc[j], U}).at(1) * h;
//   }

//   // Implicit discrete-time dynamics
//   // Function F("F", {vertcat(Xc), X0, U}, {vertcat(eq), Xf, Qf});
//   Function Feq("feq", {vertcat(Xc), X0, U}, {vertcat(eq)});
//   Function Fxf("fxf", {vertcat(Xc), X0, U}, {Xf});
//   Function Fxq("fxq", {Lc, vertcat(Xc), X0, U}, {Lc + Qf});

//   vector<double> lbg;
//   vector<double> ubg;
//   vector<double> lbx;
//   vector<double> ubx;
//   vector<double> w0;

//   // Optimization
//   MX J = 0.0;
//   vector<MX> var_xs(N + 1);
//   vector<MX> var_us(N);
//   vector<MX> var_xcs(N);

//   var_xs[0] = MX::sym("X0", 2);
//   lbx.push_back(-numeric_limits<double>::infinity());
//   lbx.push_back(-numeric_limits<double>::infinity());
//   ubx.push_back(numeric_limits<double>::infinity());
//   ubx.push_back(numeric_limits<double>::infinity());
//   w0.push_back(0);
//   w0.push_back(1);

//   Function feqmap = Feq.map(N, "openmp", 20);
//   Function fxfmap = Fxf.map(N, "openmp", 20);
//   Function fxqfold = Fxq.fold(N);

//   for (int k = 0; k < N; ++k)
//   {
//     var_us[k] = MX::sym("U_" + to_string(k));
//     lbx.push_back(-1);
//     ubx.push_back(1);
//     w0.push_back(0);

//     var_xcs[k] = MX::sym("Xc_" + to_string(k), 2 * d);
//     for (int r = 0; r < d; ++r)
//     {
//       lbx.push_back(-numeric_limits<double>::infinity());
//       lbx.push_back(-0.25);
//       ubx.push_back(numeric_limits<double>::infinity());
//       ubx.push_back(numeric_limits<double>::infinity());
//       w0.push_back(0);
//       w0.push_back(0);
//     }
//     var_xs[k + 1] = MX::sym("X_" + to_string(k + 1), 2);
//     lbx.push_back(-0.25);
//     lbx.push_back(-numeric_limits<double>::infinity());
//     ubx.push_back(numeric_limits<double>::infinity());
//     ubx.push_back(numeric_limits<double>::infinity());
//     w0.push_back(0);
//     w0.push_back(0);
//   }

//   MX xsoffset = var_xs[1];
//   MX xs = var_xs[0];
//   MX us = var_us[0];
//   MX xcs = var_xcs[0];

//   for (int k = 1; k < N; ++k)
//   {
//     xs = horzcat(xs, var_xs[k]);
//     us = horzcat(us, var_us[k]);
//     xcs = horzcat(xcs, var_xcs[k]);
//     xsoffset = horzcat(xsoffset, var_xs[k + 1]);
//   }

//   MX g = vertcat(feqmap(vector<MX>{xcs, xs, us}).at(0), fxfmap(vector<MX>{xcs, xs, us}).at(0) - xsoffset);
//   g = reshape(g, g.rows() * g.columns(), 1);
//   lbg = vector<double>(g.rows() * g.columns(), 0);
//   ubg = vector<double>(g.rows() * g.columns(), 0);
//   J = fxqfold(vector<MX>{J, xcs, xs, us}).at(0);

//   // Define the decision variables and the bounds
//   MX w = vertcat(xs, us, xcs);
//   std::cout << w << std::endl;

//   // Create the NLP solver
//   casadi::MXDict nlp = {{"x",w},
//                         {"f",0},
//                         {"g",0}};
//   Function solver = nlpsol("solver", "ipopt", nlp);
//   // solver.generate_dependencies("nlp.c");
//   // solver = nlpsol("solver", "ipopt", "nlp.c");
//   casadi::DMDict arg;
//             arg["lbg"] = lbg;
//             arg["ubg"] = ubg;
//             arg["lbx"] = lbx;
//             arg["ubx"] = ubx;
//             arg["x0"] = w0;
//   // Solve the NLP
//   DMDict sol = solver(arg);

//   return 0;
// }

// // #include <casadi/casadi.hpp>
// // #include <chrono>
// // #include <iostream>
// // #include <vector>

// // using namespace std;
// // using namespace casadi;

// // int main() {
// //   // number of inputs to evaluate in parallel
// //   int N = 300;

// //   // dummy input
// //   vector<double> dummyInput(N);
// //   for (int i = 0; i < N; ++i) {
// //     dummyInput[i] = i * (2.0 * M_PI / N);
// //   }
// //   DM dummyInputDM = DM(dummyInput);

// //   // make a dummy function that's moderately expensive to evaluate
// //   std::cout << "creating dummy function...." << endl;
// //   SX x = SX::sym("x");
// //   SX y = x;
// //   for (int k = 0; k < 100000; ++k) {
// //     y = sin(y);
// //   }
// //   Function f0 = Function("f", {x}, {y});

// //   // evaluate it serially, the old-fashioned way
// //   MX X = MX::sym("x", N);
// //   vector<MX> Y(N);
// //   for (int k = 0; k < N; ++k) {
// //     Y[k] = f0({X(k)}).at(0);
// //   }
// //   Function fNaiveParallel = Function("fParallel", {X}, {MX::vertcat(Y)});

// //   std::cout << "evaluating naive parallel function..." << endl;
// //   auto t0 = chrono::high_resolution_clock::now();
// //   vector<double> outNaive = fNaiveParallel({dummyInputDM}).at(0).nonzeros();
// //   auto t1 = chrono::high_resolution_clock::now();
// //   std::cout << "evaluated naive parallel function in " << chrono::duration<double>(t1 - t0).count() << " seconds" << endl;

// //   // evaluate it using new serial map construct
// //   Function fMap = f0.map(N);

// //   std::cout << "evaluating serial map function..." << endl;
// //   t0 = chrono::high_resolution_clock::now();
// //   vector<double> outMap = fMap({dummyInputDM}).at(0).nonzeros();
// //   t1 = chrono::high_resolution_clock::now();
// //   std::cout << "evaluated serial map function in " << chrono::duration<double>(t1 - t0).count() << " seconds" << endl;

// //   // evaluate it using new parallel map construct
// //   fMap = f0.map(N, "openmp");

// //   std::cout << "evaluating parallel map function..." << endl;
// //   t0 = chrono::high_resolution_clock::now();
// //   outMap = fMap({dummyInputDM}).at(0).nonzeros();
// //   t1 = chrono::high_resolution_clock::now();
// //   std::cout << "evaluated parallel map function in " << chrono::duration<double>(t1 - t0).count() << " seconds" << endl;

// //   return 0;
// // }

// CasADi core
#include <casadi/casadi.hpp>

using namespace casadi;

int main()
{
  // Declare variables
  SX u = SX::sym("u");                   // control
  SX r = SX::sym("r"), s = SX::sym("s"); // states
  SX x = vertcat(r, s);

  // Number of differential states
  int nx = x.size1();

  // Number of controls
  int nu = u.size1();

  // Bounds and initial guess for the control
  std::vector<double> u_min = {-0.75};
  std::vector<double> u_max = {1.0};
  std::vector<double> u_init = {0.0};

  // Bounds and initial guess for the state
  std::vector<double> x0_min = {0, 1};
  std::vector<double> x0_max = {0, 1};
  std::vector<double> x_min = {-inf, -inf};
  std::vector<double> x_max = {inf, inf};
  std::vector<double> xf_min = {0, 0};
  std::vector<double> xf_max = {0, 0};
  std::vector<double> x_init = {0, 0};

  // Final time
  double tf = 20.0;

  // Number of shooting nodes
  int ns = 100;
  double dt = tf/ns;

  // ODE right hand side and quadrature
  SX ode = vertcat((1 - s * s) * r - s + u, r);
  SX quad = r * r + s * s + u * u;
  SXDict dae = {{"x", x}, {"p", u}, {"ode", ode}, {"quad", quad}};
  // Define the system dynamics
  Function f = Function("f", {x, u}, {ode});

  // Define the Runge-Kutta 4 integration
  auto k1 = f(SXVector{x, u}).at(0);
  SX k2 = f(SXVector{x+dt/2*k1, u}).at(0);
  SX k3 = f(SXVector{x+dt/2*k2, u}).at(0);
  SX k4 = f(SXVector{x+dt*k3, u}).at(0);

  SX x_next = x + dt/6*(k1+2*k2+2*k3+k4);

  // Define the quadrature function
  Function q = Function("q", {x, u}, {quad});

  // Define the Runge-Kutta 4 integration for quadrature
  auto q1 = q(SXVector{x, u}).at(0);
  SX q2 = q(SXVector{x+dt/2*k1, u}).at(0);
  SX q3 = q(SXVector{x+dt/2*k2, u}).at(0);
  SX q4 = q(SXVector{x+dt*k3, u}).at(0);

  SX q_next = quad + dt/6*(q1+2*q2+2*q3+q4);

  Function F = Function("rk_f", {x, u}, {x_next, q_next});

  // Create an integrator (CVodes)
  // Function F = integrator("integrator", "cvodes", dae, 0, tf / ns);

  // Total number of NLP variables
  int NV = nx * (ns + 1) + nu * ns;

  // Declare variable vector for the NLP
  SX V = SX::sym("V", NV);

  // NLP variable bounds and initial guess
  std::vector<double> v_min, v_max, v_init;

  // Offset in V
  int offset = 0;

  // State at each shooting node and control for each shooting interval
  std::vector<SX> X, U;
  for (int k = 0; k < ns; ++k)
  {
    // Local state
    X.push_back(V.nz(Slice(offset, offset + nx)));
    if (k == 0)
    {
      v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
      v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
    }
    else
    {
      v_min.insert(v_min.end(), x_min.begin(), x_min.end());
      v_max.insert(v_max.end(), x_max.begin(), x_max.end());
    }
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());
    offset += nx;

    // Local control
    U.push_back(V.nz(Slice(offset, offset + nu)));
    v_min.insert(v_min.end(), u_min.begin(), u_min.end());
    v_max.insert(v_max.end(), u_max.begin(), u_max.end());
    v_init.insert(v_init.end(), u_init.begin(), u_init.end());
    offset += nu;
  }

  // State at end
  X.push_back(V.nz(Slice(offset, offset + nx)));
  v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
  v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
  v_init.insert(v_init.end(), x_init.begin(), x_init.end());
  offset += nx;

  // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
  casadi_assert(offset == NV, "");

  // Objective function
  SX J = 0;

  // Constraint function and bounds
  std::vector<SX> g;

  // Loop over shooting nodes
  for (int k = 0; k < ns; ++k)
  {
    // Create an evaluation node
    auto I_out = F(SXVector{X[k], U[k]});

    // Save continuity constraints
    g.push_back(I_out.at(0) - X[k + 1]);

    // Add objective function contribution
    J += I_out.at(1);
  }

  // NLP
  SXDict nlp = {{"x", V}, {"f", J}, {"g", vertcat(g)}};

  // Set options
  Dict opts;
  opts["ipopt.tol"] = 1e-5;
  opts["ipopt.max_iter"] = 100;
  opts["jit"] = true;
  opts["compiler"] = "shell";
  opts["jit_options.flags"] = "-fPIC -shared -O3";
  opts["jit_options.compiler"] = "gcc";

  // Create an NLP solver and buffers
  Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);

  std::map<std::string, DM> arg, res;

  // Bounds and initial guess
  arg["lbx"] = v_min;
  arg["ubx"] = v_max;
  arg["lbg"] = 0;
  arg["ubg"] = 0;
  arg["x0"] = v_init;

  // Solve the problem
  res = solver(arg);

  // Optimal solution of the NLP
  std::vector<double> V_opt(res.at("x"));

  // Get the optimal state trajectory
  std::vector<double> r_opt(ns + 1), s_opt(ns + 1);
  for (int i = 0; i <= ns; ++i)
  {
    r_opt[i] = V_opt.at(i * (nx + 1));
    s_opt[i] = V_opt.at(1 + i * (nx + 1));
  }
  // std::cout << "r_opt = " << std::endl
  //           << r_opt << std::endl;
  // std::cout << "s_opt = " << std::endl
  //           << s_opt << std::endl;

  // // Get the optimal control
  // std::vector<double> u_opt(ns);
  // for (int i = 0; i < ns; ++i)
  // {
  //   u_opt[i] = V_opt.at(nx + i * (nx + 1));
  // }
  // std::cout << "u_opt = " << std::endl
  //           << u_opt << std::endl;
}