#pragma once

#include "galileo/opt/Segment.h"
#include "galileo/opt/PseudospectralSegment.h"
#include <chrono>
#include <Eigen/Dense>

using namespace casadi;
using namespace std;

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
             * @param builders Constraint builders used to build the constraints
             * @param opts_ Options to pass to the solver
             */
            TrajectoryOpt(shared_ptr<ProblemData> problem_, vector<shared_ptr<ConstraintBuilder<ProblemData>>> builders, Dict opts_);

            /**
             * @brief Initialize the finite elements.
             *
             * @param d The degree of the finite element polynomials
             * @param X0 The initial state to deviate from
             */
            void init_finite_elements(int d, DM X0);

            /**
             * @brief Optimize and return the solution.
             *
             * @return SXVector The solution
             */
            MXVector optimize();

            /**
             * @brief Get the times where the decision variables are evaluated.
             */
            vector<double> get_global_times() const;

            /**
             * @brief Get the solution at the specified times. Times should be monotonically increasing, but do not have to be evenly spaced.
             * The solution is interpolated between the finite elements. The times should be bounded by the times given to the optimization problem.
             *
             * @param times The times to evaluate the solution at
             * @return Eigen::MatrixXd The solution at the specified times
             */
            Eigen::MatrixXd get_solution(Eigen::VectorXd &times) const;

        private:
            /**
             * @brief A Trajectory is made up of segments of finite elements.
             *
             */
            vector<shared_ptr<Segment>> trajectory;

            /**
             * @brief The solution vector.
             *
             */
            MXVector sol;

            /**
             * @brief Problem data containing constraints problem data.
             *
             */
            shared_ptr<ProblemData> problem;

            /**
             * @brief Problem data containing the objective function and dynamics.
             *
             */
            shared_ptr<GeneralProblemData> gp_data;

            /**
             * @brief Constraint builders used to build the constraints.
             *
             */
            vector<shared_ptr<ConstraintBuilder<ProblemData>>> builders;

            /**
             * @brief Casadi solver options.
             *
             */
            Dict opts;

            /**
             * @brief Nonlinear function solver.
             *
             */
            Function solver;

            /**
             * @brief Slicer to get the states.
             *
             */
            shared_ptr<States> state_indices;

            /**
             * @brief Vector of all decision variables.
             *
             */
            MXVector w;

            /**
             * @brief Vector of all constraint expressions.
             *
             */
            MXVector g;

            /**
             * @brief Vector of all general constraint lower bounds.
             *
             */
            vector<double> lbg;

            /**
             * @brief Vector of all general constraint upper bounds.
             *
             */
            vector<double> ubg;

            /**
             * @brief  Lower bounds associated with the decision variables.
             *
             */
            vector<double> lbw;

            /**
             * @brief Upper bounds associated with the decision variables.
             *
             */
            vector<double> ubw;

            /**
             * @brief Initial guess for the decision variables.
             */
            vector<double> w0;

            /**
             * @brief Expression for objective cost.
             *
             */
            MX J;

            /**
             * @brief Vector of all times where decision variables are evaluated.
             *
             */
            shared_ptr<DM> global_times;
        };

        template <class ProblemData>
        TrajectoryOpt<ProblemData>::TrajectoryOpt(shared_ptr<ProblemData> problem_, vector<shared_ptr<ConstraintBuilder<ProblemData>>> builders_, Dict opts_)
        {
            this->problem = problem_;
            this->gp_data = problem_->gp_data;
            this->builders = builders_;
            this->state_indices = problem_->states;
            this->opts = opts_;
        }

        template <class ProblemData>
        void TrajectoryOpt<ProblemData>::init_finite_elements(int d, DM X0)
        {
            assert(X0.size1() == state_indices->nx && X0.size2() == 1 && "Initial state must be a column vector");
            trajectory.clear();
            w.clear();
            g.clear();
            lbg.clear();
            ubg.clear();
            lbw.clear();
            ubw.clear();
            w0.clear();
            J = 0;

            global_times = nullptr;

            MX prev_final_state = X0;
            MX prev_final_state_deviant;
            MX curr_initial_state_deviant;

            shared_ptr<Segment> segment;

            /*Each phase should have a vector of constraint data*/
            size_t num_phases = 1;

            /*DUMMY DATA FOR TESTING*/
            vector<double> equality_back(state_indices->nx, 0.0);

            // Cheat vars just for testing the constraint framework
            SX x = SX::sym("x", state_indices->nx);
            SX u = SX::sym("u", state_indices->nu);
            SX t = SX::sym("t");

            vector<shared_ptr<ConstraintData>> G;
            shared_ptr<ConstraintData> u_bound_constraint = make_shared<ConstraintData>();

            // u_bound_constraint->global = true;
            // u_bound_constraint->upper_bound = Function("u_ubound", {t}, {1.0});
            // u_bound_constraint->lower_bound = Function("u_lbound", {t}, {-1.0});
            // u_bound_constraint->G = Function("u_bound", {x, u}, {u});
            // G.push_back(u_bound_constraint);

            /*Validation of time varying bounds*/
            shared_ptr<DecisionData> Wx = make_shared<DecisionData>();
            // Wx->upper_bound = Function("x_ubound", {t}, {SX::ones(state_indices->ndx, 1) * inf});
            // Wx->lower_bound = Function("x_lbound", {t}, {SX::ones(state_indices->ndx, 1) * -inf});
            // // Wx->lower_bound = Function("x_lbound", {t}, {SX::vertcat({-0.07 * (t - 1.0) * (t - 1.0) - 0.25, -1.0})});
            // Wx->initial_guess = Function("x_guess", {t}, {SX::zeros(state_indices->nx, 1)});
            // Wx->w = x;

            shared_ptr<DecisionData> Wu = make_shared<DecisionData>();
            // Wu->upper_bound = Function("u_ubound", {t}, {SX::ones(state_indices->nu, 1)});
            // Wu->lower_bound = Function("u_lbound", {t}, {-SX::ones(state_indices->nu, 1)});
            // Wu->initial_guess = Function("u_guess", {t}, {SX::zeros(state_indices->nu, 1)});
            // Wu->w = u;

            /*END OF DUMMY DATA*/

            printf("Starting initialization\n");
            auto start_time = chrono::high_resolution_clock::now();
            Function Phi = gp_data->Phi;

            for (size_t i = 0; i < num_phases; ++i)
            {
                /*TODO; Replace this ugly constructor with ProblemData. Most of this info should be stored in there anyways*/
                segment = make_shared<PseudospectralSegment>(gp_data, state_indices, d, 20, 1. / 20);
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
                    J += Phi(MXVector{prev_final_state}).at(0);
                }
            }
            printf("Finished initialization.\n");
            auto end_time = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = end_time - start_time;
            cout << "Elapsed time for initialization: " << elapsed.count() << endl;
        }

        template <class ProblemData>
        MXVector TrajectoryOpt<ProblemData>::optimize()
        {

            MXDict nlp = {{"x", vertcat(w)},
                                  {"f", J},
                                  {"g", vertcat(g)}};
            solver = nlpsol("solver", "ipopt", nlp, opts);

            double time_from_funcs = 0.0;
            double time_just_ipopt = 0.0;
            DMDict arg;
            arg["lbg"] = lbg;
            arg["ubg"] = ubg;
            arg["lbx"] = lbw;
            arg["ubx"] = ubw;
            arg["x0"] = w0;
            DMDict sol = solver(arg);
            w0 = sol["x"].get_elements();
            Dict stats = solver.stats();
            time_from_funcs += (double)stats["t_wall_nlp_jac_g"] + (double)stats["t_wall_nlp_hess_l"] + (double)stats["t_wall_nlp_grad_f"] + (double)stats["t_wall_nlp_g"] + (double)stats["t_wall_nlp_f"];
            time_just_ipopt += (double)stats["t_wall_total"] - time_from_funcs;
            auto tmp = MX(sol["x"]);
            cout << "Extracting solution..." << endl;

            cout << "Total seconds from Casadi functions: " << time_from_funcs << endl;
            cout << "Total seconds from Ipopt w/o function: " << time_just_ipopt << endl;
            sol = trajectory[0]->extract_solution(tmp);
            return sol;
        }

        template <class ProblemData>
        Eigen::MatrixXd TrajectoryOpt<ProblemData>::get_solution(Eigen::VectorXd &times) const
        {
            /*TODO: Only implemented for pseudospectral segments (and not very well), use visitor pattern if we want to use other types of segments*/
            int i = 0;
            auto solx = sol[0];
            int d = 0;
            Eigen::MatrixXd result(state_indices->nx, times.size());

            for (shared_ptr<Segment> seg : trajectory)
            {
                shared_ptr<PseudospectralSegment> pseg = static_pointer_cast<PseudospectralSegment>(seg);
                auto l_times = pseg->get_local_times();
                double start_time = l_times(0);
                double end_time = l_times(l_times.size1() - 1);

                /*Get the segment of times that is greater than start time and less than or equal to end time*/
                auto start_it = lower_bound(times.data(), times.data() + times.size(), start_time);
                auto end_it = upper_bound(times.data(), times.data() + times.size(), end_time);
                Eigen::VectorXd segment_times(end_it - start_it);
                copy(start_it, end_it, segment_times.data());
                DM solx_segment = MX::evalf(solx(Slice(0, solx.rows()), Slice(d, d + pseg->dX_poly.d)));
                DMVector solx_segment_vec;
                for (int i = 0; i < solx_segment.size2(); ++i)
                {
                    solx_segment_vec.push_back(solx_segment(Slice(0, solx_segment.rows()), i));
                }

                d += pseg->dX_poly.d;
                for (int j = 0; j < segment_times.size(); ++j)
                {
                    auto tmp = pseg->dX_poly.lagrange_interpolation(segment_times[j], solx_segment_vec);
                    for (int k = 0; k < tmp.rows(); ++k)
                    {
                        result(k, i) = double(tmp(k));
                    }
                    ++i;
                }
            }
            return result;
        }

        template <class ProblemData>
        vector<double> TrajectoryOpt<ProblemData>::get_global_times() const
        {
            return (*global_times).get_elements();
        }

    }
}