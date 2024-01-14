#pragma once

#include "galileo/opt/Segment.h"
#include "third-party/gnuplot-iostream/gnuplot-iostream.h"
#include <chrono>

namespace galileo
{
    namespace opt
    {
        /**
         * @brief The trajectory optimization class. This class is responsible for
            initializing the finite elements, and optimizing the trajectory.
         *
         */
        template <class ProblemData, class SegmentType>
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
            TrajectoryOpt(std::shared_ptr<ProblemData> problem_, std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders, casadi::Dict opts_);

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
             * @return casadi::SXVector The solution
             */
            casadi::MXVector optimize();

            /**
             * @brief Get the times where the decision variables are evaluated.
             */
            std::vector<double> get_global_times() const;

        private:
            /**
             * @brief A Trajectory is made up of segments of finite elements.
             *
             */
            std::vector<std::shared_ptr<SegmentType>> trajectory;

            /**
             * @brief Problem data containing the objective function and dynamics.
             *
             */
            std::shared_ptr<GeneralProblemData> gp_data;

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
            std::shared_ptr<casadi::DM> global_times;
        };

        template <class ProblemData, class SegmentType>
        TrajectoryOpt<ProblemData, SegmentType>::TrajectoryOpt(std::shared_ptr<ProblemData> problem_, std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders, casadi::Dict opts_)
        {
            this->gp_data = problem_->gp_data;
            this->state_indices = problem_->states;
            this->opts = opts_;
        }

        template <class ProblemData, class SegmentType>
        void TrajectoryOpt<ProblemData, SegmentType>::init_finite_elements(int d, casadi::DM X0)
        {
            assert(X0.size1() == this->state_indices->nx && X0.size2() == 1 && "Initial state must be a column vector");
            this->trajectory.clear();
            this->w.clear();
            this->g.clear();
            this->lbg.clear();
            this->ubg.clear();
            this->lbw.clear();
            this->ubw.clear();
            this->w0.clear();
            this->J = 0;

            this->global_times = nullptr;

            casadi::MX prev_final_state = X0;
            casadi::MX prev_final_state_deviant;
            casadi::MX curr_initial_state_deviant;

            std::shared_ptr<SegmentType> segment;

            /*Each phase should have a vector of constraint data*/
            std::size_t num_phases = 1;

            /*DUMMY DATA FOR TESTING*/
            std::vector<double> equality_back(this->state_indices->nx, 0.0);

            // Cheat vars just for testing the constraint framework
            casadi::SX x = casadi::SX::sym("x", this->state_indices->nx);
            casadi::SX u = casadi::SX::sym("u", this->state_indices->nu);
            casadi::SX t = casadi::SX::sym("t");

            std::vector<std::shared_ptr<ConstraintData>> G;
            std::shared_ptr<ConstraintData> u_bound_constraint = std::make_shared<ConstraintData>();

            // u_bound_constraint->global = true;
            // u_bound_constraint->upper_bound = casadi::Function("u_ubound", {t}, {1.0});
            // u_bound_constraint->lower_bound = casadi::Function("u_lbound", {t}, {-1.0});
            // u_bound_constraint->G = casadi::Function("u_bound", {x, u}, {u});
            // G.push_back(u_bound_constraint);

            /*Validation of time varying bounds*/
            std::shared_ptr<DecisionData> Wx = std::make_shared<DecisionData>();
            // Wx->upper_bound = casadi::Function("x_ubound", {t}, {casadi::SX::ones(this->state_indices->ndx, 1) * inf});
            // Wx->lower_bound = casadi::Function("x_lbound", {t}, {casadi::SX::ones(this->state_indices->ndx, 1) * -inf});
            // // Wx->lower_bound = casadi::Function("x_lbound", {t}, {casadi::SX::vertcat({-0.07 * (t - 1.0) * (t - 1.0) - 0.25, -1.0})});
            // Wx->initial_guess = casadi::Function("x_guess", {t}, {casadi::SX::zeros(this->state_indices->nx, 1)});
            // Wx->w = x;

            std::shared_ptr<DecisionData> Wu = std::make_shared<DecisionData>();
            // Wu->upper_bound = casadi::Function("u_ubound", {t}, {casadi::SX::ones(this->state_indices->nu, 1)});
            // Wu->lower_bound = casadi::Function("u_lbound", {t}, {-casadi::SX::ones(this->state_indices->nu, 1)});
            // Wu->initial_guess = casadi::Function("u_guess", {t}, {casadi::SX::zeros(this->state_indices->nu, 1)});
            // Wu->w = u;

            /*END OF DUMMY DATA*/

            printf("Starting initialization\n");
            auto start_time = std::chrono::high_resolution_clock::now();
            Function Phi = this->gp_data->Phi;

            for (std::size_t i = 0; i < num_phases; ++i)
            {
                /*TODO; Replace this ugly constructor with ProblemData. Most of this info should be stored in there anyways*/
                segment = std::make_shared<SegmentType>(this->gp_data, d, 20, 1. / 20, this->global_times, this->state_indices, X0, prev_final_state, G, Wx, Wu, this->J, this->w, this->g);

                this->trajectory.push_back(segment);

                segment->fill_lbg_ubg(this->lbg, this->ubg);
                segment->fill_lbw_ubw(this->lbw, this->ubw);
                segment->fill_w0(this->w0);
                this->global_times = segment->get_global_times();

                /*Initial state constraint*/
                if (i == 0)
                {
                    auto curr_initial_state = segment->get_initial_state();
                    this->g.push_back(prev_final_state - curr_initial_state);
                    this->lbg.insert(this->lbg.end(), equality_back.begin(), equality_back.end());
                    this->ubg.insert(this->ubg.end(), equality_back.begin(), equality_back.end());
                }
                /*Continuity constraint for the state deviant between phases*/
                else if (i > 0)
                {
                    curr_initial_state_deviant = segment->get_initial_state_deviant();
                    /*For general jump map functions you can use the following syntax:*/
                    // g.push_back(jump_map_function(casadi::MXVector{prev_final_state_deviant, curr_initial_state_deviant}).at(0));
                    this->g.push_back(prev_final_state_deviant - curr_initial_state_deviant);
                    this->lbg.insert(this->lbg.end(), equality_back.begin(), equality_back.end() - 1);
                    this->ubg.insert(this->ubg.end(), equality_back.begin(), equality_back.end() - 1);
                }
                prev_final_state = segment->get_final_state();
                prev_final_state_deviant = segment->get_final_state_deviant();

                /*Terminal cost*/
                if (i == num_phases - 1 && i != 0)
                {
                    this->J += Phi(casadi::MXVector{prev_final_state}).at(0);
                }
            }
            printf("Finished initialization.\n");
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            std::cout << "Elapsed time for initialization: " << elapsed.count() << std::endl;
        }

        template <class ProblemData, class SegmentType>
        casadi::MXVector TrajectoryOpt<ProblemData, SegmentType>::optimize()
        {
            // std::cout << "w: " << vertcat(this->w) << std::endl;
            // std::cout << "g: " << vertcat(this->g) << std::endl;
            // std::cout << "size w: " << vertcat(this->w).size() << std::endl;
            // std::cout << "size g: " << vertcat(this->g).size() << std::endl;
            // std::cout << this->lb.size() << std::endl;
            // std::cout << this->ub.size() << std::endl;

            casadi::MXDict nlp = {{"x", vertcat(this->w)},
                                  {"f", this->J},
                                  {"g", vertcat(this->g)}};
            this->solver = casadi::nlpsol("solver", "ipopt", nlp, this->opts);
            // auto grad_f = casadi::Function("grad_f", {vertcat(this->w)}, {casadi::MX::gradient(this->J, vertcat(this->w))}).expand();
            // auto jac_g = casadi::Function("jac_g", {vertcat(this->w)}, {casadi::MX::jacobian(vertcat(this->g), vertcat(this->w))}).expand();
            // auto eval_grad_f = grad_f(vertcat(this->w));
            // auto gradeval = eval_grad_f.at(0);
            // std::cout << "eval_grad_f: " << gradeval(1215,0) << std::endl;

            // auto eval_jac_g = jac_g(vertcat(this->w));
            // auto jaceval = eval_jac_g.at(0);
            // std::cout << "eval_jac_g: " << jaceval << std::endl;

            double time_from_funcs = 0.0;
            double time_just_ipopt = 0.0;
            casadi::DMDict arg;
            arg["lbg"] = this->lbg;
            arg["ubg"] = this->ubg;
            arg["lbx"] = this->lbw;
            arg["ubx"] = this->ubw;
            arg["x0"] = this->w0;
            casadi::DMDict sol = this->solver(arg);
            this->w0 = sol["x"].get_elements();
            casadi::Dict stats = this->solver.stats();
            time_from_funcs += (double)stats["t_wall_nlp_jac_g"] + (double)stats["t_wall_nlp_hess_l"] + (double)stats["t_wall_nlp_grad_f"] + (double)stats["t_wall_nlp_g"] + (double)stats["t_wall_nlp_f"];
            time_just_ipopt += (double)stats["t_wall_total"] - time_from_funcs;
            auto tmp = casadi::MX(sol["x"]);
            std::cout << "Extracting solution..." << std::endl;

            std::cout << "Total seconds from Casadi functions: " << time_from_funcs << std::endl;
            std::cout << "Total seconds from Ipopt w/o function: " << time_just_ipopt << std::endl;
            return this->trajectory[0]->extract_solution(tmp);
        }

        template <class ProblemData, class SegmentType>
        std::vector<double> TrajectoryOpt<ProblemData, SegmentType>::get_global_times() const
        {
            return (*this->global_times).get_elements();
        }

    }
}