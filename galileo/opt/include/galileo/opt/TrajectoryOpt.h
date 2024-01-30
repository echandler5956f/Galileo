#pragma once

#include "galileo/opt/Segment.h"
#include "galileo/opt/PseudospectralSegment.h"
#include <chrono>
#include <Eigen/Dense>

namespace galileo
{
    namespace opt
    {
        /**
         * @brief The trajectory optimization class. This class is responsible for
            initializing the finite elements, and optimizing the trajectory.
         *
         */
        template <class ProblemData>
        class TrajectoryOpt
        {
        public:
            /**
             * @brief Construct a new Trajectory Opt object.
             *
             * @param problem_ Problem data containing the objective function and dynamics
             * @param builders_ Constraint builders used to build the constraints
             * @param opts_ Options to pass to the solver
             */
            TrajectoryOpt(std::shared_ptr<ProblemData> problem_, std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders_, casadi::Dict opts_);

            /**
             * @brief Initialize the finite elements.
             *
             * @param d The degree of the finite element polynomials
             * @param X0 The initial state to deviate from
             */
            void init_finite_elements(int d, casadi::DM X0);

            /**
             * @brief Optimize and return the solution.
             *
             * @return SXVector The solution
             */
            casadi::MXVector optimize();

            /**
             * @brief Get the times where the decision variables are evaluated.
             */
            std::vector<double> get_global_times() const;

            /**
             * @brief Get the solution at the specified times. Times should be monotonically increasing, but do not have to be evenly spaced.
             * The solution is interpolated between the finite elements. The times should be bounded by the times given to the optimization problem.
             *
             * @param times The times to evaluate the solution at
             * @return Eigen::MatrixXd The solution at the specified times
             */
            Eigen::MatrixXd get_solution(Eigen::VectorXd &times) const;

            Eigen::VectorXd get_segment_times(Eigen::VectorXd &times, double initial_time, double end_time) const;

            void process_segment_times(std::vector<double> &l_times_vec, Eigen::VectorXd &segment_times, casadi::DM &solx_segment, int degree, const std::shared_ptr<LagrangePolynomial> poly, Eigen::MatrixXd &result, int i, int j) const;

        private:
            /**
             * @brief A Trajectory is made up of segments of finite elements.
             *
             */
            std::vector<std::shared_ptr<Segment>> trajectory;

            /**
             * @brief The solution vector.
             *
             */
            casadi::MXVector sol;

            /**
             * @brief Problem data containing constraints problem data.
             *
             */
            std::shared_ptr<ProblemData> problem;

            /**
             * @brief Problem data containing the objective function and dynamics.
             *
             */
            std::shared_ptr<GeneralProblemData> gp_data;

            /**
             * @brief Constraint builders used to build the constraints.
             *
             */
            std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders;

            /**
             * @brief Casadi solver options.
             *
             */
            casadi::Dict opts;

            /**
             * @brief Nonlinear function solver.
             *
             */
            casadi::Function solver;

            /**
             * @brief Slicer to get the states.
             *
             */
            std::shared_ptr<States> state_indices;

            /**
             * @brief Vector of all decision variables.
             *
             */
            casadi::MXVector w;

            /**
             * @brief Vector of all constraint expressions.
             *
             */
            casadi::MXVector g;

            /**
             * @brief Vector of all general constraint lower bounds.
             *
             */
            std::vector<double> lbg;

            /**
             * @brief Vector of all general constraint upper bounds.
             *
             */
            std::vector<double> ubg;

            /**
             * @brief  Lower bounds associated with the decision variables.
             *
             */
            std::vector<double> lbw;

            /**
             * @brief Upper bounds associated with the decision variables.
             *
             */
            std::vector<double> ubw;

            /**
             * @brief Initial guess for the decision variables.
             */
            std::vector<double> w0;

            /**
             * @brief Expression for objective cost.
             *
             */
            casadi::MX J;

            /**
             * @brief Vector of all times where decision variables are evaluated.
             *
             */
            casadi::DM global_times;
        };

        template <class ProblemData>
        TrajectoryOpt<ProblemData>::TrajectoryOpt(std::shared_ptr<ProblemData> problem_, std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders_, casadi::Dict opts_)
        {
            this->problem = problem_;
            this->gp_data = problem_->gp_data;
            this->builders = builders_;
            this->state_indices = problem_->states;
            this->opts = opts_;
        }

        template <class ProblemData>
        void TrajectoryOpt<ProblemData>::init_finite_elements(int d, casadi::DM X0)
        {
            assert(X0.size1() == state_indices->nx && X0.size2() == 1 && "Initial state must be a column vector");
            trajectory.clear();
            global_times = casadi::DM(0, 0);
            w.clear();
            g.clear();
            lbg.clear();
            ubg.clear();
            lbw.clear();
            ubw.clear();
            w0.clear();
            J = 0;

            casadi::MX prev_final_state = X0;
            casadi::MX prev_final_state_deviant;
            casadi::MX curr_initial_state_deviant;

            std::shared_ptr<Segment> segment;

            /*Each phase should have a vector of constraint data*/
            size_t num_phases = 1;

            /*DUMMY DATA FOR TESTING*/
            std::vector<double> equality_back(state_indices->nx, 0.0);

            // Cheat vars just for testing the constraint framework
            casadi::SX x = casadi::SX::sym("x", state_indices->nx);
            casadi::SX u = casadi::SX::sym("u", state_indices->nu);
            casadi::SX t = casadi::SX::sym("t");

            std::vector<std::shared_ptr<ConstraintData>> G;
            std::shared_ptr<ConstraintData> u_bound_constraint = std::make_shared<ConstraintData>();

            // u_bound_constraint->global = true;
            // u_bound_constraint->upper_bound = Function("u_ubound", {t}, {1.0});
            // u_bound_constraint->lower_bound = Function("u_lbound", {t}, {-1.0});
            // u_bound_constraint->G = Function("u_bound", {x, u}, {u});
            // G.push_back(u_bound_constraint);

            /*Validation of time varying bounds*/
            std::shared_ptr<DecisionData> Wx = std::make_shared<DecisionData>();
            // Wx->upper_bound = Function("x_ubound", {t}, {SX::ones(state_indices->ndx, 1) * inf});
            // Wx->lower_bound = Function("x_lbound", {t}, {SX::ones(state_indices->ndx, 1) * -inf});
            // // Wx->lower_bound = Function("x_lbound", {t}, {SX::vertcat({-0.07 * (t - 1.0) * (t - 1.0) - 0.25, -1.0})});
            // Wx->initial_guess = Function("x_guess", {t}, {SX::zeros(state_indices->nx, 1)});
            // Wx->w = x;

            std::shared_ptr<DecisionData> Wu = std::make_shared<DecisionData>();
            // Wu->upper_bound = Function("u_ubound", {t}, {SX::ones(state_indices->nu, 1)});
            // Wu->lower_bound = Function("u_lbound", {t}, {-SX::ones(state_indices->nu, 1)});
            // Wu->initial_guess = Function("u_guess", {t}, {SX::zeros(state_indices->nu, 1)});
            // Wu->w = u;

            /*END OF DUMMY DATA*/

            printf("Starting initialization\n");
            auto start_time = std::chrono::high_resolution_clock::now();
            casadi::Function Phi = gp_data->Phi;

            for (size_t i = 0; i < num_phases; ++i)
            {
                /*TODO; Replace this ugly constructor with ProblemData. Most of this info should be stored in there anyways*/
                segment = std::make_shared<PseudospectralSegment>(gp_data, state_indices, d, 20, 10. / 20);
                segment->initialize_local_time_vector(global_times);
                segment->initialize_knot_segments(X0, prev_final_state);
                segment->initialize_expression_graph(G, Wx, Wu);
                segment->evaluate_expression_graph(J, w, g);

                trajectory.push_back(segment);

                segment->fill_lbg_ubg(lbg, ubg);
                segment->fill_lbw_ubw(lbw, ubw);
                segment->fill_w0(w0);

                /*Initial state constraint*/
                if (i == 0)
                {
                    auto curr_initial_state = segment->get_initial_state();
                    g.push_back(prev_final_state - curr_initial_state);
                    lbg.insert(lbg.end(), equality_back.begin(), equality_back.end());
                    ubg.insert(ubg.end(), equality_back.begin(), equality_back.end());
                }
                /*Continuity constraint for the state deviant between phases*/
                else if (i > 0)
                {
                    curr_initial_state_deviant = segment->get_initial_state_deviant();
                    /*For general jump map functions you can use the following syntax:*/
                    // g.push_back(jump_map_function(MXVector{prev_final_state_deviant, curr_initial_state_deviant}).at(0));
                    g.push_back(prev_final_state_deviant - curr_initial_state_deviant);
                    lbg.insert(lbg.end(), equality_back.begin(), equality_back.end() - 1);
                    ubg.insert(ubg.end(), equality_back.begin(), equality_back.end() - 1);
                }
                prev_final_state = segment->get_final_state();
                prev_final_state_deviant = segment->get_final_state_deviant();

                /*Terminal cost*/
                if (i == num_phases - 1)
                {
                    J += Phi(casadi::MXVector{prev_final_state}).at(0);
                }
            }
            printf("Finished initialization.\n");
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            std::cout << "Elapsed time for initialization: " << elapsed.count() << std::endl;
        }

        template <class ProblemData>
        casadi::MXVector TrajectoryOpt<ProblemData>::optimize()
        {

            casadi::MXDict nlp = {{"x", vertcat(w)},
                                  {"f", J},
                                  {"g", vertcat(g)}};
            solver = casadi::nlpsol("solver", "ipopt", nlp, opts);

            double time_from_funcs = 0.0;
            double time_just_ipopt = 0.0;
            casadi::DMDict arg;
            arg["lbg"] = lbg;
            arg["ubg"] = ubg;
            arg["lbx"] = lbw;
            arg["ubx"] = ubw;
            arg["x0"] = w0;
            casadi::DMDict result = solver(arg);
            w0 = result["x"].get_elements();
            casadi::Dict stats = solver.stats();
            time_from_funcs += (double)stats["t_wall_nlp_jac_g"] + (double)stats["t_wall_nlp_hess_l"] + (double)stats["t_wall_nlp_grad_f"] + (double)stats["t_wall_nlp_g"] + (double)stats["t_wall_nlp_f"];
            time_just_ipopt += (double)stats["t_wall_total"] - time_from_funcs;
            auto tmp = casadi::MX(result["x"]);
            std::cout << "Extracting solution..." << std::endl;

            std::cout << "Total seconds from Casadi functions: " << time_from_funcs << std::endl;
            std::cout << "Total seconds from Ipopt w/o function: " << time_just_ipopt << std::endl;
            sol = trajectory[0]->extract_solution(tmp);
            return sol;
        }

        template <class ProblemData>
        Eigen::MatrixXd TrajectoryOpt<ProblemData>::get_solution(Eigen::VectorXd &times) const
        {
            auto start_time = std::chrono::high_resolution_clock::now();
            int i = 0;
            auto solx = sol[0];
            int count = 0;
            Eigen::MatrixXd result(state_indices->nx, times.size());

            for (std::shared_ptr<Segment> seg : trajectory)
            {
                std::shared_ptr<PseudospectralSegment> pseg = std::dynamic_pointer_cast<PseudospectralSegment>(seg);
                std::vector<double> l_times_vec = pseg->get_local_times().get_elements();
                double initial_time = l_times_vec[0];
                double end_time = l_times_vec[l_times_vec.size() - 1];
                int degree = pseg->get_degree() + 1;
                int num_knots = pseg->get_knot_num();
                auto polynomial = pseg->get_dXPoly();

                Eigen::VectorXd segment_times = get_segment_times(times, initial_time, end_time);
                casadi::DM solx_segment = casadi::MX::evalf(solx(casadi::Slice(0, solx.rows()), casadi::Slice(count, count + degree * num_knots)));

                for (Eigen::Index j = 0; j < segment_times.size(); ++j)
                {
                    process_segment_times(l_times_vec, segment_times, solx_segment, degree, polynomial, result, i, j);
                    ++i;
                }
                count += degree * num_knots;
            }
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            std::cout << "Elapsed time for get_solution: " << elapsed.count() << std::endl;
            return result;
        }

        template <class ProblemData>
        Eigen::VectorXd TrajectoryOpt<ProblemData>::get_segment_times(Eigen::VectorXd &times, double initial_time, double end_time) const
        {
            auto start_it = std::lower_bound(times.data(), times.data() + times.size(), initial_time);
            auto end_it = std::upper_bound(times.data(), times.data() + times.size(), end_time);
            Eigen::VectorXd segment_times(end_it - start_it);
            std::copy(start_it, end_it, segment_times.data());
            return segment_times;
        }

        template <class ProblemData>
        void TrajectoryOpt<ProblemData>::process_segment_times(std::vector<double> &l_times_vec, Eigen::VectorXd &segment_times, casadi::DM &solx_segment, int degree, const std::shared_ptr<LagrangePolynomial> poly, Eigen::MatrixXd &result, int i, int j) const
        {
            auto it = std::lower_bound(l_times_vec.data(), l_times_vec.data() + l_times_vec.size(), segment_times(j));
            int index = (it - l_times_vec.data() - 1) / (degree);
            casadi::DM solx_knot_segment = solx_segment(casadi::Slice(0, solx_segment.rows()), casadi::Slice(index * degree, (index * degree) + degree));
            casadi::DMVector solx_segment_vec;
            for (casadi_int k = 0; k < solx_knot_segment.size2(); ++k)
            {
                solx_segment_vec.push_back(solx_knot_segment(casadi::Slice(0, solx_knot_segment.rows()), k));
            }
            double scaled_time = (segment_times[j] - l_times_vec[index * degree]) / (l_times_vec[(index * degree) + degree] - l_times_vec[index * degree]);
            std::vector<double> tmp = poly->barycentric_interpolation(scaled_time, solx_segment_vec).get_elements();
            for (std::size_t k = 0; k < tmp.size(); ++k)
            {
                result(k, i) = tmp[k];
            }
        }

        template <class ProblemData>
        std::vector<double> TrajectoryOpt<ProblemData>::get_global_times() const
        {
            return global_times.get_elements();
        }
    }
}