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

            // this->w0 = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., X0(0,0).get_elements().at(0), X0(1,0).get_elements().at(0), 0., 0., 0., 0., 0., 0., 0., 0., 0.};

            std::vector<double> equality_back(this->state_indices->nx, 0.0);
            std::shared_ptr<PseudospectralSegment> ps;
            std::vector<std::shared_ptr<ConstraintData>> G;
            // std::size_t i = 0;
            printf("Starting\n");
            auto startinit = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < 1; ++i)
            {
                auto start = std::chrono::high_resolution_clock::now();
                ps = std::make_shared<PseudospectralSegment>(d, 1000, 10. / 1000, this->state_indices, this->Fint);
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                std::cout << "Time taken for ps Constructor call: " << duration.count() << " microseconds" << std::endl;

                start = std::chrono::high_resolution_clock::now();
                ps->initialize_knot_segments(prev_final_state);
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                std::cout << "Time taken for initialize_knot_segments call: " << duration.count() << " microseconds" << std::endl;

                /*TODO: Fill with user defined functions, and handle global/phase-dependent/time-varying constraints*/
                start = std::chrono::high_resolution_clock::now();
                ps->initialize_expression_graph(this->F, this->L, G);
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                std::cout << "Time taken for initialize_expression_graph call: " << duration.count() << " microseconds" << std::endl;

                this->trajectory.push_back(ps);
                start = std::chrono::high_resolution_clock::now();
                ps->evaluate_expression_graph(this->J, this->w, this->g);
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                std::cout << "Time taken for evaluate_expression_graph call: " << duration.count() << " microseconds" << std::endl;

                start = std::chrono::high_resolution_clock::now();
                ps->fill_lbg_ubg(this->lbg, this->ubg);
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                std::cout << "Time taken for fill_lbg_ubg call: " << duration.count() << " microseconds" << std::endl;

                start = std::chrono::high_resolution_clock::now();
                ps->fill_times(this->all_times);
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                std::cout << "Time taken for fill_times call: " << duration.count() << " microseconds" << std::endl;

                if (i == 0)
                {
                    auto curr_initial_state = ps->get_initial_state();
                    this->g.push_back(prev_final_state - curr_initial_state);
                    this->lbg.insert(this->lbg.end(), equality_back.begin(), equality_back.end());
                    this->ubg.insert(this->ubg.end(), equality_back.begin(), equality_back.end());
                }
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
                if (i == 2)
                {
                    /*Add terminal cost*/
                    this->J += this->Phi(casadi::SXVector{prev_final_state}).at(0);
                }
                ++i;
            }
            auto endinit = std::chrono::high_resolution_clock::now();
            auto durationinit = std::chrono::duration_cast<std::chrono::microseconds>(endinit - startinit);
            std::cout << "Time taken for all Pseudospectral segment methods: " << durationinit.count() << " microseconds" << std::endl;
            printf("Finished init\n");
        }

        casadi::SXVector TrajectoryOpt::optimize()
        {
            std::cout << "w: " << vertcat(this->w) << std::endl;
            // std::cout << "g: " << vertcat(this->g) << std::endl;
            // std::cout << "size w: " << vertcat(this->w).size() << std::endl;
            // std::cout << "size g: " << vertcat(this->g).size() << std::endl;
            // std::cout << this->lb.size() << std::endl;
            // std::cout << this->ub.size() << std::endl;
            // printf("HERE!!!!!!!!!!!!!!!!!!\n");

            casadi::SXDict nlp = {{"x", vertcat(this->w)},
                                  {"f", this->J},
                                  {"g", vertcat(this->g)}};
            this->solver = casadi::nlpsol("solver", "ipopt", nlp, this->opts);
            casadi::DMDict arg;
            arg["lbg"] = this->lbg;
            arg["ubg"] = this->ubg;
            // arg["lbx"] = this->lbw;
            // arg["ubx"] = this->ubw;
            arg["x0"] = this->w0;
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