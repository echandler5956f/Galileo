#include "galileo/opt/TrajectoryOpt.h"

namespace galileo
{
    namespace opt
    {
        TrajectoryOpt::TrajectoryOpt(casadi::Dict opts_, std::shared_ptr<States> state_indices_, std::shared_ptr<GeneralProblemData> problem)
        {
            this->opts = opts_;
            this->state_indices = state_indices_;

            this->Fint = problem->Fint;
            this->Fdif = problem->Fdif;
            this->F = problem->F;
            this->L = problem->L;
            this->Phi = problem->Phi;
        }

        void TrajectoryOpt::init_finite_elements(int d, casadi::DM X0)
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

            std::shared_ptr<PseudospectralSegment> ps;

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
            // Wx->upper_bound = casadi::Function("x_ubound", {t}, {casadi::SX::ones(this->state_indices->ndx, 1)});
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
            for (std::size_t i = 0; i < num_phases; ++i)
            {
                ps = std::make_shared<PseudospectralSegment>(d, 20, 1. / 20, this->global_times, this->state_indices, this->Fint, this->Fdif);
                ps->initialize_knot_segments(X0, prev_final_state);

                /*TODO: Fill with user defined functions, and handle global/phase-dependent/time-varying constraints*/
                /*TODO: Pass in W which holds the decision variable bounds and initial guess data*/
                ps->initialize_expression_graph(this->F, this->L, G, Wx, Wu);
                std::cout << "Finished initializing expression graph" << std::endl;

                this->trajectory.push_back(ps);
                ps->evaluate_expression_graph(this->J, this->w, this->g);

                ps->fill_lbg_ubg(this->lbg, this->ubg);
                ps->fill_lbw_ubw(this->lbw, this->ubw);
                ps->fill_w0(this->w0);
                this->global_times = ps->get_global_times();

                /*Initial state constraint*/
                if (i == 0)
                {
                    auto curr_initial_state = ps->get_initial_state();
                    this->g.push_back(prev_final_state - curr_initial_state);
                    this->lbg.insert(this->lbg.end(), equality_back.begin(), equality_back.end());
                    this->ubg.insert(this->ubg.end(), equality_back.begin(), equality_back.end());
                }
                /*Continuity constraint for the state deviant between phases*/
                else if (i > 0)
                {
                    curr_initial_state_deviant = ps->get_initial_state_deviant();
                    /*For general jump map functions you can use the following syntax:*/
                    // g.push_back(jump_map_function(casadi::MXVector{prev_final_state_deviant, curr_initial_state_deviant}).at(0));
                    this->g.push_back(prev_final_state_deviant - curr_initial_state_deviant);
                    this->lbg.insert(this->lbg.end(), equality_back.begin(), equality_back.end() - 1);
                    this->ubg.insert(this->ubg.end(), equality_back.begin(), equality_back.end() - 1);
                }
                prev_final_state = ps->get_final_state();
                prev_final_state_deviant = ps->get_final_state_deviant();

                /*Terminal cost*/
                if (i == num_phases - 1 && i != 0)
                {
                    this->J += this->Phi(casadi::MXVector{prev_final_state}).at(0);
                }
            }
            printf("Finished initialization.\n");
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            std::cout << "Elapsed time for initialization: " << elapsed.count() << std::endl;
        }

        casadi::MXVector TrajectoryOpt::optimize()
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

        std::vector<double> TrajectoryOpt::get_global_times() const
        {
            return (*this->global_times).get_elements();
        }
    }
}