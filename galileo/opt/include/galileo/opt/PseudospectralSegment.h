#pragma once

#include "galileo/opt/Segment.h"

namespace galileo
{
    namespace opt
    {

        /**
         * @brief Helper class for storing polynomial information.
         *
         */
        class LagrangePolynomial
        {
        public:
            /**
             * @brief Construct a new Lagrange Polynomial object.
             *
             */
            LagrangePolynomial(){};

            /**
             * @brief Construct a new Lagrange Polynomial object. Compute and store the coefficients for a given degree and collocation scheme.
             *
             * @param d_ Degree of the polynomial
             * @param scheme Collocation scheme: "radau" or "legendre"
             */
            LagrangePolynomial(int d_, const std::string &scheme = "radau");

            /**
             * @brief Perform symbolic Lagrange Interpolation, which, given a time from the Lagrange time scale, interpolates terms to find the value at time t.
             *
             * @param t Time to interpolate at
             * @param terms Terms at knot points to use for interpolation
             * @return const casadi::SX Resultant expression for the symbolic interpolated value
             */
            template <typename Scalar>
            Scalar lagrange_interpolation(double t, const std::vector<Scalar> terms) const;

            casadi::DM barycentric_interpolation(double t, const casadi::DMVector terms) const;

            /**
             * @brief Degree of the polynomial.
             *
             */
            int d;

            /**
             * @brief The roots of the polynomial.
             *
             */
            std::vector<double> tau_root;

            /**
             * @brief Quadrature coefficients.
             *
             */
            std::vector<double> B;

            /**
             * @brief Collocation coeffficients.
             *
             */
            std::vector<std::vector<double>> C;

            /**
             * @brief Continuity coefficients.
             *
             */
            std::vector<double> D;
        };

        /**
         * @brief PseudospectalSegment class.
         *
         */
        class PseudospectralSegment : public Segment
        {
        public:
            /**
             * @brief Construct a new Pseudospectral Segment object.
             *
             * @param problem Pointer to the problem data
             * @param st_m_ Pointer to the state indices helper
             * @param d Polynomial degree
             * @param knot_num_ Number of knots in the segment
             * @param h_ Period of each knot segment
             *
             */
            PseudospectralSegment(std::shared_ptr<GeneralProblemData> problem, std::shared_ptr<States> st_m_, int d, int knot_num_, double h_);

            /**
             * @brief Initialize the relevant expressions.
             *
             * @param d Polynomial degree
             */
            void initialize_expression_variables(int d);

            /**
             * @brief Initialize the vector of local times which constraints are evaluated at.
             *
             * @param global_times Vector of global times
             */
            void initialize_local_time_vector(casadi::DM &global_times);

            /**
             * @brief Initialize the vector of times which coincide to the decision variables U occur at.
             *
             */
            void initialize_u_time_vector();

            /**
             * @brief Create all the knot segments.
             *
             * @param x0_global Global starting state to integrate from (used for initial guess)
             * @param x0_local Local starting state to integrate from
             *
             */
            void initialize_knot_segments(casadi::DM x0_global, casadi::MX x0_local);

            /**
             * @brief Build the function graph.
             *
             * @param G Vector of constraint data
             * @param Wx Decision bound and initial guess data for the state
             * @param Wu Decision bound and initial guess data for the input
             */
            void initialize_expression_graph(std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu);

            /**
             * @brief Evaluate the expressions with the actual decision variables.
             *
             * @param J0 Accumulated cost so far
             * @param w Decision variable vector to fill
             * @param g Constraint vector to fill
             */
            void evaluate_expression_graph(casadi::MX &J0, casadi::MXVector &w, casadi::MXVector &g);

            /**
             * @brief Extract the solution from the decision variable vector.
             *
             * @param w Decision variable vector
             * @return casadi::MX Solution values
             */
            casadi::MXVector extract_solution(casadi::MX &w) const;

            /**
             * @brief Get the initial state.
             *
             * @return casadi::MX The initial state
             */
            casadi::MX get_initial_state() const;

            /**
             * @brief Get the initial state deviant.
             *
             * @return casadi::MX The initial state deviant
             */
            casadi::MX get_initial_state_deviant() const;

            /**
             * @brief Get the final state deviant.
             *
             * @return casadi::MX The final state deviant
             */
            casadi::MX get_final_state_deviant() const;

            /**
             * @brief Get the actual final state.
             *
             * @return casadi::MX The final state.
             */
            casadi::MX get_final_state() const;

            /**
             * @brief Get the local times vector.
             *
             * @return casadi::DM The local times vector
             */
            casadi::DM get_local_times() const;

            /**
             * @brief Fills the lower bounds on decision variable (lbw) and upper bounds on decision variable (ubw) vectors with values.
             *
             * This function takes in two vectors, lbw and ubw, and fills them with values.
             * The filled values represent the constraint lower and upper bounds on decision variables.
             *
             * @param lbw The vector to be filled with lower bound values on decision variables.
             * @param ubw The vector to be filled with upper bound values on decision variables.
             */
            void fill_lbw_ubw(std::vector<double> &lbw, std::vector<double> &ubw);

            /**
             * @brief Fills the lower bounds on general constraints (lbg) and upper bounds on general constraints (ubg) vectors with values.
             *
             * This function takes in two vectors, lbg and ubg, and fills them with values.
             * The filled values represent the general constraint lower and upper bounds.
             *
             * @param lbg The vector to be filled with general lower bound values.
             * @param ubg The vector to be filled with general upper bound values.
             */
            void fill_lbg_ubg(std::vector<double> &lbg, std::vector<double> &ubg);

            /**
             * @brief Fills the initial guess vector (w0) with values.
             *
             * This function takes in a vector, w0, and fills it with values.
             * The filled values represent the initial guess for the decision variables.
             *
             * @param w0 The vector to be filled with initial guess values.
             */
            void fill_w0(std::vector<double> &w0) const;

            /**
             * @brief Returns the starting and ending index in w.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_decision_variables() const;

            /**
             * @brief Returns the starting and ending index in g (call after evaluate_expression_graph!).
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_constraint_expressions() const;

            /**
             * @brief Returns the starting and ending index in lbg/ubg.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_constraint_bounds() const;

            /**
             * @brief Returns the starting and ending index in lbw/ubw.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_decision_bounds() const;

            /**
             * @brief Get the knot num
             *
             * @return int
             */
            int get_knot_num() const
            {
                return knot_num;
            }

            /**
             * @brief Input polynomial. Helper object to store polynomial information for the input.
             *
             */
            LagrangePolynomial U_poly;

            /**
             * @brief State polynomial. Helper object to store polynomial information for the state.
             *
             */
            LagrangePolynomial dX_poly;

        private:
            /**
             * @brief Helper function to process a vector of type MX.
             *
             * This function takes a vector of type casadi::MX and performs some processing on it.
             * It creates a temporary vector by copying the input vector and removing the last element.
             * The modified vector is then returned.
             *
             * @param vec The input vector of type casadi::MX.
             * @return The processed vector of type casadi::MX.
             */
            casadi::MX processVector(casadi::MXVector &vec) const;

            /**
             * @brief Process the offset vector by removing the first element and concatenating the remaining elements horizontally.
             *
             * @param vec The input offset vector.
             * @return The processed offset vector.
             */
            casadi::MX processOffsetVector(casadi::MXVector &vec) const;

            /**
             * @brief Collocation state decision variables.
             *
             */
            casadi::MXVector dXc_var_vec;

            /**
             * @brief Collocation input decision variables.
             *
             */
            casadi::MXVector U_var_vec;

            /**
             * @brief Collocation input decision expressions at the state collocation points.
             *
             * Decision variables of control and state are potentially approximated by different degree polynomials.
             *
             */
            casadi::SXVector U_at_c_vec;

            /**
             * @brief Collocation state decision expressions at the collocation points.
             *
             */
            casadi::SXVector x_at_c_vec;

            /**
             * @brief Knot point deviants state decision variables.
             *
             */
            casadi::MXVector dX0_var_vec;

            /**
             * @brief Knot point state expressions (integral functions of the deviants).
             *
             */
            casadi::MXVector X0_var_vec;

            /**
             * @brief casadi::Function map for extracting the solution from the ocp solution vector.
             *
             */
            casadi::Function sol_map_func;

            /**
             * @brief Solution function.
             *
             */
            casadi::Function get_sol_func;

            /**
             * @brief Implicit discrete-time function map. This function map returns the vector of collocation equations
                necessary to match the derivative defect between the approximated dynamics and actual system
                dynamics.
             *
             */
            casadi::Function collocation_constraint_map;

            /**
             * @brief Implicit discrete-time function map. The map which matches the approximated final state expression with the initial
                state of the next segment.
             *
             */
            casadi::Function xf_constraint_map;

            /**
             * @brief Implicit discrete-time function map. The accumulated cost across all the knot segments found using quadrature rules.
             *
             */
            casadi::Function q_cost_fold;

            /**
             * @brief User defined constraints, which are functions with certain bounds associated with them.
             *
             */
            std::vector<casadi::Function> general_constraint_maps;

            /**
             * @brief Lower bounds associated with the general constraint maps.
             *
             */
            casadi::DM general_lbg;

            /**
             * @brief Upper bounds associated with the general constraint maps.
             *
             */
            casadi::DM general_ubg;

            /**
             * @brief Lower bounds associated with the decision variable constraint maps.
             *
             */
            casadi::DM general_lbw;

            /**
             * @brief Upper bounds associated with the decision variable constraint maps.
             *
             */
            casadi::DM general_ubw;

            /**
             * @brief Initial guess for associated with this segment
             */
            casadi::DM w0;

            /**
             * @brief Integrator function.
             *
             */
            casadi::Function Fint;

            /**
             * @brief Difference function.
             *
             */
            casadi::Function Fdif;

            /**
             * @brief Dynamics function.
             *
             */
            casadi::Function F;

            /**
             * @brief Cost function.
             *
             */
            casadi::Function L;

            /**
             * @brief Collocation states used to build the expression graphs.
             *
             */
            casadi::SXVector dXc;

            /**
             * @brief Collocation inputs used to build the expression graphs.
             *
             */
            casadi::SXVector Uc;

            /**
             * @brief Knot states deviants used to build the expression graphs.
             *
             */
            casadi::SX dX0;

            /**
             * @brief Knot states used to build the expression graphs.
             *
             */
            casadi::SX X0;

            /**
             * @brief Accumulator expression used to build the expression graphs.
             *
             */
            casadi::SX Lc;

            /**
             * @brief Helper for indexing the state variables.
             *
             */
            std::shared_ptr<States> st_m;

            /**
             * @brief Number of knot segments.
             *
             */
            int knot_num;

            /**
             * @brief Vector of local times including knot points. Note that this coincides with the times for the decision variables of x.
             *
             */
            casadi::DM local_times;

            /**
             * @brief Vector of unique times for this segment, including terminal but not including initial.
             * Used for bound constraint evaluation, so that we don't overconstrain variables which occur at the same time
             * e.g, x0 and xf of adjacent segments
             *
             * In the frame of the global times (NOT the local times)
             *
             */
            casadi::DM collocation_times;

            /**
             * @brief Vector of times for the knot points for this segment.
             *
             * In the frame of the global times (NOT the local times)
             */
            casadi::DM knot_times;

            /**
             * @brief Vector of times for the decision variables of u.
             */
            casadi::DM u_times;
        };

        LagrangePolynomial::LagrangePolynomial(int d_, const std::string &scheme)
        {
            d = d_;
            /*Choose collocation points*/
            tau_root = casadi::collocation_points(d, scheme);
            tau_root.insert(tau_root.begin(), 0);

            /*Coefficients of the quadrature function*/
            B.resize(d + 1);

            /*Coefficients of the collocation equation*/
            C.resize(d + 1);
            for (int j = 0; j < d + 1; ++j)
                C[j].resize(d + 1);

            /*Coefficients of the continuity equation*/
            D.resize(d + 1);

            /*For all collocation points*/
            for (int j = 0; j < d + 1; ++j)
            {
                /*Construct Lagrange polynomials to get the polynomial basis at the collocation point*/
                casadi::Polynomial p = 1;
                for (int r = 0; r < d + 1; ++r)
                {
                    if (r != j)
                    {
                        p *= casadi::Polynomial(-tau_root[r], 1) / (tau_root[j] - tau_root[r]);
                    }
                }
                /*Evaluate the polynomial at the final time to get the coefficients of the continuity equation*/
                D[j] = p(1.0);

                /*Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation*/
                casadi::Polynomial dp = p.derivative();
                for (int r = 0; r < d + 1; ++r)
                {
                    C[j][r] = dp(tau_root[r]);
                }
                casadi::Polynomial pint = p.anti_derivative();
                B[j] = pint(1.0);
            }
        }

        template <typename Scalar>
        Scalar LagrangePolynomial::lagrange_interpolation(double t, const std::vector<Scalar> terms) const
        {
            assert((t >= 0.0) && (t <= 1.0) && "t must be in the range [0,1]");
            Scalar result = 0;
            for (int j = 0; j < d; ++j)
            {
                Scalar term = terms[j];
                for (int r = 0; r < d + 1; ++r)
                {
                    if (r != j)
                    {
                        term *= (t - tau_root[r]) / (tau_root[j] - tau_root[r]);
                    }
                }
                result += term;
            }
            return result;
        }

        casadi::DM LagrangePolynomial::barycentric_interpolation(double t, const casadi::DMVector terms) const
        {
            casadi::DMVector w;
            for (int j = 0; j < tau_root.size(); ++j)
            {
                w.push_back(1.0);
            }

            // Compute the barycentric weights
            for (int j = 0; j < d; ++j)
            {
                for (int r = 0; r < tau_root.size(); ++r)
                {
                    if (r != j)
                    {
                        w[j] *= (tau_root[j] - tau_root[r]);
                    }
                }
                w[j] = 1.0 / w[j];
            }

            // Compute the interpolated value
            casadi::DM numerator = 0.0, denominator = 0.0;
            for (int i = 0; i < tau_root.size(); ++i)
            {
                if (std::abs(t - tau_root[i]) < 1e-6)
                {
                    return terms[i];
                }
                casadi::DM term = w[i] / (t - tau_root[i]);
                numerator += term * terms[i];
                denominator += term;
            }

            return numerator / denominator;
        }

        PseudospectralSegment::PseudospectralSegment(std::shared_ptr<GeneralProblemData> problem, std::shared_ptr<States> st_m_, int d, int knot_num_, double h_)
        {
            auto Fint_ = problem->Fint;
            auto Fdif_ = problem->Fdif;
            auto F_ = problem->F;
            auto L_ = problem->L;

            assert(d > 0 && d < 10 && "d must be greater than 0 and less than 10");
            assert(h_ > 0 && "h must be a positive duration");
            assert(Fint_.n_in() == 3 && "Fint must have 3 inputs");
            assert(Fint_.n_out() == 1 && "Fint must have 1 output");
            assert(Fdif_.n_in() == 3 && "Fdif must have 3 inputs");
            assert(Fdif_.n_out() == 1 && "Fdif must have 1 output");

            Fint_.assert_size_in(0, st_m_->nx, 1);
            Fint_.assert_size_in(1, st_m_->ndx, 1);
            Fint_.assert_size_in(2, 1, 1);
            Fint_.assert_size_out(0, st_m_->nx, 1);

            Fdif_.assert_size_in(0, st_m_->nx, 1);
            Fdif_.assert_size_in(1, st_m_->nx, 1);
            Fdif_.assert_size_in(2, 1, 1);
            Fdif_.assert_size_out(0, st_m_->ndx, 1);

            this->knot_num = knot_num_;
            this->h = h_;
            this->st_m = st_m_;
            this->Fint = Fint_;
            this->Fdif = Fdif_;
            this->F = F_;
            this->L = L_;
            this->T = (knot_num)*h;

            initialize_expression_variables(d);
        }

        void PseudospectralSegment::initialize_expression_variables(int d)
        {
            dXc.clear();
            Uc.clear();

            dX_poly = LagrangePolynomial(d);
            U_poly = LagrangePolynomial(1);

            for (int j = 0; j < dX_poly.d; ++j)
            {
                dXc.push_back(casadi::SX::sym("dXc_" + std::to_string(j), st_m->ndx, 1));
                if (j < U_poly.d)
                {
                    Uc.push_back(casadi::SX::sym("Uc_" + std::to_string(j), st_m->nu, 1));
                }
            }
            dX0 = casadi::SX::sym("dX0", st_m->ndx, 1);
            X0 = casadi::SX::sym("X0", st_m->nx, 1);
            Lc = casadi::SX::sym("Lc", 1, 1);
        }

        void PseudospectralSegment::initialize_local_time_vector(casadi::DM &global_times)
        {
            local_times = casadi::DM::zeros(knot_num * (dX_poly.d + 1), 1);
            collocation_times = casadi::DM::zeros(knot_num * (dX_poly.d), 1);
            knot_times = casadi::DM::zeros(knot_num + 1, 1);
            int i = 0;
            int unique_i = 0;
            int j = 0;
            double kh = 0.0;
            for (int k = 0; k < knot_num; ++k)
            {
                kh = k * h;
                local_times(i, 0) = kh;
                knot_times(k, 0) = kh;
                ++i;
                for (j = 0; j < dX_poly.d; ++j)
                {
                    local_times(i, 0) = kh + dX_poly.tau_root[j + 1] * h;
                    if (i > 0 && unique_i < collocation_times.size1())
                    {
                        collocation_times(unique_i, 0) = local_times(i, 0);
                        ++unique_i;
                    }
                    ++i;
                }
            }
            knot_times(knot_times.size1() - 2, 0) = T - h;
            knot_times(knot_times.size1() - 1, 0) = T;

            double start_time = 0.0;
            if (global_times.is_empty() == false)
            {
                start_time = global_times(global_times.size1() - 1, 0).get_elements()[0];
                local_times += start_time;
                global_times = vertcat(global_times, local_times);
            }
            else
            {
                global_times = local_times;
            }
            collocation_times += start_time;
            knot_times += start_time;
            initialize_u_time_vector();
        }

        void PseudospectralSegment::initialize_u_time_vector()
        {
            u_times = casadi::DM::zeros(knot_num * (U_poly.d), 1);
            int j = 0;
            double kh = 0.0;
            u_times(0, 0) = kh;
            int i = 1;
            for (int k = 0; k < knot_num; ++k)
            {
                kh = k * h;
                for (j = 0; j < U_poly.d; ++j)
                {
                    if (i < u_times.size1())
                        u_times(i, 0) = kh + U_poly.tau_root[j + 1] * h;
                    else
                        break;
                    ++i;
                }
            }
        }

        void PseudospectralSegment::initialize_knot_segments(casadi::DM x0_global_, casadi::MX x0_local_)
        {
            x0_global = x0_global_;
            x0_local = x0_local_;
            assert(x0_local.size1() == st_m->nx && x0_local.size2() == 1 && "x0 must be a column std::vector of size nx");

            dXc_var_vec.clear();
            U_var_vec.clear();
            dX0_var_vec.clear();
            X0_var_vec.clear();
            U_at_c_vec.clear();
            x_at_c_vec.clear();
            for (int k = 0; k < knot_num; ++k)
            {
                dXc_var_vec.push_back(casadi::MX::sym("dXc_" + std::to_string(k), st_m->ndx * dX_poly.d, 1));
                U_var_vec.push_back(casadi::MX::sym("U_" + std::to_string(k), st_m->nu * U_poly.d, 1));
            }

            for (int k = 0; k < knot_num + 1; ++k)
            {
                dX0_var_vec.push_back(casadi::MX::sym("dX0_" + std::to_string(k), st_m->ndx, 1));
                X0_var_vec.push_back(Fint(casadi::MXVector{x0_local, dX0_var_vec[k], 1.0}).at(0));
            }
        }

        void PseudospectralSegment::initialize_expression_graph(std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu)
        {
            assert(F.n_in() == 2 && "F must have 2 inputs");
            assert(F.n_out() == 1 && "F must have 1 output");

            assert(L.n_in() == 2 && "L must have 2 inputs");
            assert(L.n_out() == 1 && "L must have 1 output");

            F.assert_size_in(0, st_m->nx, 1);
            F.assert_size_in(1, st_m->nu, 1);
            F.assert_size_out(0, st_m->ndx, 1);

            L.assert_size_in(0, st_m->nx, 1);
            L.assert_size_in(1, st_m->nu, 1);
            L.assert_size_out(0, 1, 1);

            /*Collocation equations*/
            casadi::SXVector eq;
            /*State at the end of the collocation interval*/
            casadi::SX dXf = dX_poly.D[0] * dX0;
            /*Cost at the end of the collocation interval*/
            casadi::SX Qf = 0;
            /*Actual state at collocation points*/
            casadi::SXVector x_at_c;
            /*U interpolated at the dx polynomial collocation points*/
            casadi::SXVector u_at_c;
            casadi::SXVector tmp_x;
            casadi::SXVector tmp_dx;
            tmp_x.push_back(X0);
            tmp_dx.push_back(dX0);

            for (int j = 0; j < dX_poly.d; ++j)
            {
                double dt_j = (dX_poly.tau_root[j + 1] - dX_poly.tau_root[j]) * h;
                /*Expression for the state derivative at the collocation point*/
                casadi::SX dxp = dX_poly.C[0][j + 1] * dX0;
                for (int r = 0; r < dX_poly.d; ++r)
                {
                    dxp += dX_poly.C[r + 1][j + 1] * dXc[r];
                }
                /*dXc must exist in a Euclidean space, but we need x_c in order to evaluate the objective. Fint can simply return dXc[j] if the states are already Euclidean*/
                casadi::SX x_c = Fint(casadi::SXVector{X0, dXc[j], dt_j}).at(0);
                // casadi::SX u_c = U_poly.lagrange_interpolation(dX_poly.tau_root[j], Uc);
                casadi::SX u_c = Uc[0];

                x_at_c.push_back(x_c);
                u_at_c.push_back(u_c);
                tmp_x.push_back(x_c);
                tmp_dx.push_back(dXc[j]);

                /*Append collocation equations*/
                eq.push_back(h * F(casadi::SXVector{x_c, u_c}).at(0) - dxp);

                /*Add cost contribution*/
                casadi::SXVector L_out = L(casadi::SXVector{x_c, u_c});
                /*This is fine as long as the cost is not related to the Lie Group elements. See the state integrator and dX for clarity*/
                Qf += dX_poly.B[j + 1] * L_out.at(0) * h;
                // Qf += U_poly.B(j + 1) * L_out.at(1) * h;

                dXf += dX_poly.D[j + 1] * dXc[j];
            }

            casadi::Dict opts;
            // // opts["cse"] = true;
            // opts["jit"] = true;
            // opts["jit_options.flags"] = "-Ofast -march=native -ffast-math";
            // opts["jit_options.compiler"] = "gcc";
            // // opts["jit_options.temp_suffix"] = false;
            // opts["compiler"] = "shell";

            auto collocation_constraint = casadi::Function("feq",
                                                           casadi::SXVector{X0, vertcat(dXc), dX0, vertcat(Uc)},
                                                           casadi::SXVector{vertcat(eq)}, opts);
            // collocation_constraint.generate("feq");
            // int flag1 = system("gcc -fPIC -shared -O3 feq.c -o feq.so");
            // casadi_assert(flag1==0, "Compilation failed");
            // auto collocation_constraint_new = external("feq");

            // auto collocation_constraint_adj1 = collocation_constraint.reverse(1);
            // collocation_constraint_adj1.generate("adj1_feq");
            // int flag2 = system("gcc -fPIC -shared -O3 adj1_feq.c -o adj1_feq.so");
            // casadi_assert(flag2==0, "Compilation failed");
            // auto collocation_constraint_adj = external("adj1_feq");

            auto xf_constraint = casadi::Function("fxf",
                                                  casadi::SXVector{X0, vertcat(dXc), dX0, vertcat(Uc)},
                                                  casadi::SXVector{dXf}, opts);
            // xf_constraint.generate("fxf");
            // int flag2 = system("gcc -fPIC -shared -O3 fxf.c -o fxf.so");
            // casadi_assert(flag2==0, "Compilation failed");
            // xf_constraint = external("fxf");

            auto q_cost = casadi::Function("fxq", casadi::SXVector{Lc, X0, vertcat(dXc), dX0, vertcat(Uc)},
                                           casadi::SXVector{Lc + Qf}, opts);
            // q_cost.generate("fxq");
            // int flag3 = system("gcc -fPIC -shared -O3 fxq.c -o fxq.so");
            // casadi_assert(flag3==0, "Compilation failed");
            // q_cost = external("fxq");

            /*Implicit discrete-time equations*/
            collocation_constraint_map = collocation_constraint.map(knot_num, "openmp");
            /*When you evaluate this map, subtract by the knot points list offset by 1 to be correct*/
            xf_constraint_map = xf_constraint.map(knot_num, "openmp");
            q_cost_fold = q_cost.fold(knot_num);

            sol_map_func = casadi::Function("sol_map",
                                            casadi::SXVector{X0, vertcat(dXc), dX0, vertcat(Uc)},
                                            casadi::SXVector{horzcat(tmp_x), horzcat(u_at_c)})
                               .map(knot_num, "serial");

            casadi_int N = collocation_constraint_map.size1_out(0) * collocation_constraint_map.size2_out(0) +
                           xf_constraint_map.size1_out(0) * xf_constraint_map.size2_out(0);
            auto tmp = N;

            std::vector<tuple_size_t> ranges_G;

            /*Map the constraint to each collocation point, and then map the mapped constraint to each knot segment*/
            for (size_t i = 0; i < G.size(); ++i)
            {
                auto g_data = G[i];

                assert(g_data->G.n_in() == 2 && "G must have 2 inputs");
                g_data->G.assert_size_in(0, st_m->nx, 1);
                g_data->G.assert_size_in(1, st_m->nu, 1);
                /*TODO: Add assertions to check the bounds functions here!!!*/

                auto tmap = casadi::Function("fg",
                                             casadi::SXVector{X0, vertcat(dXc), dX0, vertcat(Uc)},
                                             casadi::SXVector{vertcat(g_data->G.map(dX_poly.d, "serial")((casadi::SXVector{horzcat(x_at_c), horzcat(u_at_c)})))})
                                .map(knot_num, "serial");
                general_constraint_maps.push_back(tmap);
                ranges_G.push_back(tuple_size_t(N, N + tmap.size1_out(0) * tmap.size2_out(0)));
                N += tmap.size1_out(0) * tmap.size2_out(0);
            }

            general_lbg.resize(N, 1);
            general_ubg.resize(N, 1);
            general_lbg(casadi::Slice(0, tmp)) = casadi::DM::zeros(tmp, 1);
            general_ubg(casadi::Slice(0, tmp)) = casadi::DM::zeros(tmp, 1);

            for (casadi_int i = 0; i < G.size(); ++i)
            {
                auto g_data = G[i];
                general_lbg(casadi::Slice(std::get<0>(ranges_G[i])), std::get<1>(ranges_G[i])) =
                    vertcat(g_data->lower_bound.map(knot_num * (dX_poly.d), "serial")(collocation_times));
                general_ubg(casadi::Slice(std::get<0>(ranges_G[i])), std::get<1>(ranges_G[i])) =
                    vertcat(g_data->upper_bound.map(knot_num * (dX_poly.d), "serial")(collocation_times));
            }

            auto Ndxknot = st_m->ndx * (knot_num + 1);
            auto Ndx = st_m->ndx * (dX_poly.d + 1) * knot_num + st_m->ndx;
            auto Ndxcol = Ndx - Ndxknot;

            auto Nu = st_m->nu * U_poly.d * knot_num;
            w0 = casadi::DM::zeros(Ndx + Nu, 1);
            general_lbw = -casadi::DM::inf(Ndx + Nu, 1);
            general_ubw = casadi::DM::inf(Ndx + Nu, 1);

            /*Transform initial guess for x to an initial guess for dx, using f_dif, the inverse of f_int*/

            casadi::MX xkg_sym = casadi::MX::sym("xkg", st_m->nx, 1);
            casadi::MX xckg_sym = casadi::MX::sym("Xckg", st_m->nx * dX_poly.d, 1);

            if (!Wx->initial_guess.is_null())
            {
                auto xg = vertcat(Wx->initial_guess.map(knot_num + 1, "serial")(knot_times));
                casadi::Function dxg_func = casadi::Function("xg_fun", casadi::MXVector{xkg_sym}, casadi::MXVector{Fdif(casadi::MXVector{x0_global, xkg_sym, 1.0}).at(0)})
                                                .map(knot_num + 1, "serial");
                w0(casadi::Slice(0, Ndxknot)) = reshape(dxg_func(casadi::DMVector{xg}).at(0), Ndxknot, 1);
                /*The transformation of xc to dxc is a slightly less trivial. While x_k = fint(x0_init, dx_k), for xc_k, we have xc_k = fint(x_k, dxc_k) which is equivalent to xc_k = fint(fint(x0_init, dx_k), dxc_k).
                Thus, dxc_k = fdif(fint(x0_init, dx_k), xc_k)). This could be done with maps like above, but it is not necessary.*/
                auto xc_g = vertcat(Wx->initial_guess.map((dX_poly.d) * knot_num, "serial")(collocation_times));
                for (casadi_int i = 0; i < knot_num; ++i)
                {
                    auto xk = xg(casadi::Slice(i * st_m->nx, (i + 1) * st_m->nx));
                    auto xck = xc_g(casadi::Slice(i * st_m->nx * dX_poly.d, (i + 1) * st_m->nx * dX_poly.d));
                    for (casadi_int j = 0; j < dX_poly.d; ++j)
                    {
                        w0(casadi::Slice(Ndxknot + i * st_m->ndx * dX_poly.d + j * st_m->ndx, Ndxknot + i * st_m->ndx * dX_poly.d + (j + 1) * st_m->ndx)) =
                            reshape(Fdif(casadi::DMVector{xk, xck(casadi::Slice(j * st_m->nx, (j + 1) * st_m->nx)), h}).at(0), st_m->ndx, 1);
                    }
                }
            }

            if (!Wx->lower_bound.is_null() && !Wx->upper_bound.is_null())
            {
                general_lbw(casadi::Slice(0, Ndxknot)) = reshape(vertcat(Wx->lower_bound.map(knot_num + 1, "serial")(knot_times)), Ndxknot, 1);
                general_ubw(casadi::Slice(0, Ndxknot)) = reshape(vertcat(Wx->upper_bound.map(knot_num + 1, "serial")(knot_times)), Ndxknot, 1);
                general_lbw(casadi::Slice(Ndxknot, Ndx)) = reshape(vertcat(Wx->lower_bound.map((dX_poly.d) * knot_num, "serial")(collocation_times)), Ndxcol, 1);
                general_ubw(casadi::Slice(Ndxknot, Ndx)) = reshape(vertcat(Wx->upper_bound.map((dX_poly.d) * knot_num, "serial")(collocation_times)), Ndxcol, 1);
            }

            if (!Wu->initial_guess.is_null())
            {
                w0(casadi::Slice(Ndx, Ndx + Nu)) = reshape(vertcat(Wu->initial_guess.map(U_poly.d * knot_num, "serial")(u_times)), Nu, 1);
            }

            if (!Wu->lower_bound.is_null() && !Wu->upper_bound.is_null())
            {
                general_lbw(casadi::Slice(Ndx, Ndx + Nu)) = reshape(vertcat(Wu->lower_bound.map(U_poly.d * knot_num, "serial")(u_times)), Nu, 1);
                general_ubw(casadi::Slice(Ndx, Ndx + Nu)) = reshape(vertcat(Wu->upper_bound.map(U_poly.d * knot_num, "serial")(u_times)), Nu, 1);
            }
        }

        casadi::MX PseudospectralSegment::processVector(casadi::MXVector &vec) const
        {
            casadi::MXVector temp = vec;
            temp.pop_back();
            return horzcat(temp);
        }

        casadi::MX PseudospectralSegment::processOffsetVector(casadi::MXVector &vec) const
        {
            casadi::MXVector temp = vec;
            temp.erase(temp.begin());
            return horzcat(temp);
        }

        void PseudospectralSegment::evaluate_expression_graph(casadi::MX &J0, casadi::MXVector &w, casadi::MXVector &g)
        {
            assert(J0.size1() == 1 && J0.size2() == 1 && "J0 must be a scalar");

            casadi::MXVector result;

            casadi::MX xs = processVector(X0_var_vec);
            casadi::MX dxs = processVector(dX0_var_vec);
            casadi::MX dxcs = horzcat(dXc_var_vec);
            casadi::MX us = horzcat(U_var_vec);
            casadi::MX xs_offset = processOffsetVector(X0_var_vec);
            casadi::MX dxs_offset = processOffsetVector(dX0_var_vec);

            casadi::MXVector solmap_restult = sol_map_func(casadi::MXVector{xs, dxcs, dxs, us});
            casadi::MX all_xs = solmap_restult.at(0);
            casadi::MX all_us = solmap_restult.at(1);

            /*This section cannot get much faster, it is bounded by the time to evaluate the constraint*/
            casadi::MX col_con_mat = collocation_constraint_map(casadi::MXVector{xs, dxcs, dxs, us}).at(0);
            casadi::MX xf_con_mat = xf_constraint_map(casadi::MXVector{xs, dxcs, dxs, us}).at(0);
            dxs_offset = reshape(dxs_offset, dxs_offset.size1() * dxs_offset.size2(), 1);

            result.push_back(reshape(col_con_mat, col_con_mat.size1() * col_con_mat.size2(), 1));
            result.push_back(reshape(xf_con_mat, xf_con_mat.size1() * xf_con_mat.size2(), 1) -
                             dxs_offset);

            for (size_t i = 0; i < general_constraint_maps.size(); ++i)
            {
                casadi::MX g_con_mat = general_constraint_maps[i](casadi::MXVector{xs, dxcs, dxs, us}).at(0);
                result.push_back(reshape(g_con_mat, g_con_mat.size1() * g_con_mat.size2(), 1));
            }

            casadi::MX cost = q_cost_fold(casadi::MXVector{J0, xs, dxcs, dxs, us}).at(0);
            J0 = cost;
            /*where g of this segment starts*/
            size_t g_size = g.size();
            /*Use move to avoid copying the vectors. Reserve space for g in advance outside of PseudospectralSegment.*/
            g.insert(g.end(), make_move_iterator(result.begin()), make_move_iterator(result.end()));
            g_range = tuple_size_t(g_size, g.size());

            /*where w of this segment starts*/
            size_t w_size = w.size();
            /*Use move to avoid copying the vectors. Reserve space for w in advance outside of PseudospectralSegment.*/
            w.insert(w.end(), make_move_iterator(dX0_var_vec.begin()), make_move_iterator(dX0_var_vec.end()));
            w.insert(w.end(), make_move_iterator(dXc_var_vec.begin()), make_move_iterator(dXc_var_vec.end()));
            w.insert(w.end(), make_move_iterator(U_var_vec.begin()), make_move_iterator(U_var_vec.end()));

            w_range = tuple_size_t(w_size, accumulate(w.begin() + w_size, w.end(), 0.0, [](int sum, const casadi::MX &item)
                                                      { return sum + item.size1() * item.size2(); }));
            get_sol_func = casadi::Function("func",
                                            casadi::MXVector({vertcat(casadi::MXVector(w.begin() + w_size, w.begin() + w.size()))}),
                                            casadi::MXVector({all_xs, all_us}));
        }

        casadi::MXVector PseudospectralSegment::extract_solution(casadi::MX &w) const
        {
            return get_sol_func(casadi::MXVector{w(casadi::Slice(casadi_int(std::get<0>(w_range)), casadi_int(std::get<1>(w_range))))});
        }

        casadi::MX PseudospectralSegment::get_initial_state_deviant() const
        {
            return dX0_var_vec.front();
        }

        casadi::MX PseudospectralSegment::get_initial_state() const
        {
            return X0_var_vec.front();
        }

        casadi::MX PseudospectralSegment::get_final_state_deviant() const
        {
            return dX0_var_vec.back();
        }

        casadi::MX PseudospectralSegment::get_final_state() const
        {
            return X0_var_vec.back();
        }

        casadi::DM PseudospectralSegment::get_local_times() const
        {
            return local_times;
        }

        void PseudospectralSegment::fill_lbw_ubw(std::vector<double> &lbw, std::vector<double> &ubw)
        {
            /*where lb/ub of this segment starts*/
            auto bw_size = lbw.size();
            std::vector<double> element_access1 = general_lbw.get_elements();
            std::vector<double> element_access2 = general_ubw.get_elements();

            lbw.insert(lbw.end(), element_access1.begin(), element_access1.end());
            ubw.insert(ubw.end(), element_access2.begin(), element_access2.end());
            lbw_ubw_range = tuple_size_t(bw_size, lbw.size());
        }

        void PseudospectralSegment::fill_lbg_ubg(std::vector<double> &lbg, std::vector<double> &ubg)
        {
            /*where lb/ub of this segment starts*/
            auto bg_size = lbg.size();
            std::vector<double> element_access1 = general_lbg.get_elements();
            std::vector<double> element_access2 = general_ubg.get_elements();

            lbg.insert(lbg.end(), element_access1.begin(), element_access1.end());
            ubg.insert(ubg.end(), element_access2.begin(), element_access2.end());

            lbg_ubg_range = tuple_size_t(bg_size, lbg.size());
        }

        void PseudospectralSegment::fill_w0(std::vector<double> &all_w0) const
        {
            std::vector<double> element_access1 = w0.get_elements();
            all_w0.insert(all_w0.end(), element_access1.begin(), element_access1.end());
        }

        tuple_size_t PseudospectralSegment::get_range_idx_decision_variables() const
        {
            return w_range;
        }

        tuple_size_t PseudospectralSegment::get_range_idx_constraint_expressions() const
        {
            return g_range;
        }

        tuple_size_t PseudospectralSegment::get_range_idx_constraint_bounds() const
        {
            return lbg_ubg_range;
        }

        tuple_size_t PseudospectralSegment::get_range_idx_decision_bounds() const
        {
            return lbw_ubw_range;
        }

    }
}