#include "galileo/variables/TrajectoryOpt.h"

namespace galileo
{
    namespace variables
    {
        TrajectoryOpt::TrajectoryOpt(casadi::Dict opts_, std::shared_ptr<States> state_indices_, std::shared_ptr<ProblemData> problem)
        {
            this->opts = opts_;
            this->state_indices = state_indices_;

            this->Fint = problem->Fint;
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

            this->all_times.clear();

            casadi::SX prev_final_state = X0;
            casadi::SX prev_final_state_deviant;
            casadi::SX curr_initial_state_deviant;

            std::vector<double> equality_back(this->state_indices->nx, 0.0);
            std::shared_ptr<PseudospectralSegment> ps;

            /*Each phase should have a vector of constraint data*/
            std::size_t num_phases = 1;
            std::vector<std::shared_ptr<ConstraintData>> G;
            std::shared_ptr<DecisionData> W;

            printf("Starting initialization\n");
            for (std::size_t i = 0; i < num_phases; ++i)
            {
                ps = std::make_shared<PseudospectralSegment>(d, 2, 1. / 2, this->state_indices, this->Fint);
                ps->initialize_knot_segments(prev_final_state);

                /*TODO: Fill with user defined functions, and handle global/phase-dependent/time-varying constraints*/
                /*TODO: Pass in W which holds the decision variable bounds and initial guess data*/
                ps->initialize_expression_graph(this->F, this->L, G, W);

                this->trajectory.push_back(ps);
                ps->evaluate_expression_graph(this->J, this->w, this->g);

                ps->fill_lbg_ubg(this->lbg, this->ubg);
                ps->fill_lbw_ubw(this->lbw, this->ubw);
                ps->fill_w0(this->w0);
                ps->fill_times(this->all_times);

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
                    // g.push_back(jump_map_function(casadi::SXVector{prev_final_state_deviant, curr_initial_state_deviant}).at(0));
                    this->g.push_back(prev_final_state_deviant - curr_initial_state_deviant);
                    this->lbg.insert(this->lbg.end(), equality_back.begin(), equality_back.end() - 1);
                    this->ubg.insert(this->ubg.end(), equality_back.begin(), equality_back.end() - 1);
                }
                prev_final_state = ps->get_final_state();
                prev_final_state_deviant = ps->get_final_state_deviant();

                /*Terminal cost*/
                if (i == num_phases - 1)
                {
                    this->J += this->Phi(casadi::SXVector{prev_final_state}).at(0);
                }
            }
            printf("Finished initialization.\n");
        }

        casadi::SXVector TrajectoryOpt::optimize()
        {
            std::cout << "w: " << vertcat(this->w) << std::endl;
            std::cout << "g: " << vertcat(this->g) << std::endl;
            // std::cout << "size w: " << vertcat(this->w).size() << std::endl;
            // std::cout << "size g: " << vertcat(this->g).size() << std::endl;
            // std::cout << this->lb.size() << std::endl;
            // std::cout << this->ub.size() << std::endl;

            casadi::SXDict nlp = {{"x", vertcat(this->w)},
                                  {"f", this->J},
                                  {"g", vertcat(this->g)}};
            this->solver = casadi::nlpsol("solver", "ipopt", nlp, this->opts);
            casadi::DMDict arg;
            arg["lbg"] = this->lbg;
            arg["ubg"] = this->ubg;
            // arg["lbx"] = this->lbw;
            // arg["ubx"] = this->ubw;
            // arg["x0"] = this->w0;
            auto sol = this->solver(arg);
            auto tmp = casadi::SX(sol["x"]);
            std::cout << "Extracting solution..." << std::endl;
            return this->trajectory[0]->extract_solution(tmp);
        }

        std::vector<double> TrajectoryOpt::get_all_times()
        {
            return this->all_times;
        }
    }
}