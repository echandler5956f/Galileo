#include "galileo/variables/TrajectoryOpt.h"

// namespace acro
// {
//     namespace variables
//     {
//         template<class ProblemData>
//         TrajectoryOpt<ProblemData>::TrajectoryOpt(casadi::Dict opts_, States *state_indices_, ProblemData *problem)
//         {
//             this->opts = opts_;
//             this->state_indices = state_indices_;

//             this->Fint = problem->Fint;
//             this->F = problem->F;
//             this->L = problem->L;
//             this->Phi = problem->Phi;
//         }

//         template<class ProblemData>
//         void TrajectoryOpt<ProblemData>::init_finite_elements(int d, casadi::DM X0)
//         {
//             this->w.clear();
//             this->g.clear();
//             this->lb.clear();
//             this->ub.clear();
//             this->J = 0;

//             this->all_times.clear();

//             casadi::SX prev_final_state = X0;
//             casadi::SX prev_final_state_deviant;
//             casadi::SX curr_initial_state_deviant;

//             std::vector<double> equality_back(this->state_indices->nx, 0.0);
//             // std::size_t i = 0;
//             printf("Starting\n");
//             for (std::size_t i = 0; i < 1; ++i)
//             {
//                 auto ps = PseudospectralSegment(d, 2, 0.5, this->state_indices, this->Fint);
//                 ps.initialize_knot_segments(prev_final_state);
//                 /*TODO: Fill with user defined functions, and handle global/phase-dependent/time-varying constraints*/
//                 std::vector<std::shared_ptr<ConstraintData>> G;
//                 ps.initialize_expression_graph(this->F, this->L, G);
//                 this->trajectory.push_back(ps);
//                 ps.evaluate_expression_graph(this->J, this->g);
//                 ps.fill_lb_ub(this->lb, this->ub);
//                 ps.fill_w(this->w);
//                 ps.fill_times(this->all_times);
//                 if (i == 0)
//                 {
//                     auto curr_initial_state = ps.get_initial_state();
//                     this->g.push_back(prev_final_state - curr_initial_state);
//                     this->lb.insert(this->lb.end(), equality_back.begin(), equality_back.end());
//                     this->ub.insert(this->ub.end(), equality_back.begin(), equality_back.end());
//                 }
//                 else if (i > 0)
//                 {
//                     curr_initial_state_deviant = ps.get_initial_state_deviant();
//                     /*For general jump map functions you can use the following syntax:*/
//                     // g.push_back(jump_map_function(casadi::SXVector{prev_final_state_deviant, curr_initial_state_deviant}).at(0));
//                     this->g.push_back(prev_final_state_deviant - curr_initial_state_deviant);
//                     this->lb.insert(this->lb.end(), equality_back.begin(), equality_back.end()-1);
//                     this->ub.insert(this->ub.end(), equality_back.begin(), equality_back.end()-1);
//                 }
//                 prev_final_state = ps.get_final_state();
//                 prev_final_state_deviant = ps.get_final_state_deviant();
//                 if (i == 2)
//                 {
//                     /*Add terminal cost*/
//                     this->J += this->Phi(casadi::SXVector{prev_final_state}).at(0);
//                 }
//                 ++i;
//             }
//             printf("Finished init\n");
//         }

//         template<class ProblemData>
//         casadi::DMDict TrajectoryOpt<ProblemData>::optimize()
//         {
//             std::cout << "w: " << vertcat(this->w) << std::endl;
//             std::cout << "g: " << vertcat(this->g) << std::endl;
//             std::cout << "size w: " << vertcat(this->w).size() << std::endl;
//             std::cout << "size g: " << vertcat(this->g).size() << std::endl;
//             std::cout << this->lb.size() << std::endl;
//             std::cout << this->ub.size() << std::endl;
//             printf("HERE!!!!!!!!!!!!!!!!!!\n");

//             casadi::SXDict nlp = {{"x", vertcat(this->w)},
//                                   {"f", this->J},
//                                   {"g", vertcat(this->g)}};
//             this->solver = casadi::nlpsol("solver", "ipopt", nlp, this->opts);
//             casadi::DMDict arg;
//             arg["lbg"] = this->lb;
//             arg["ubg"] = this->ub;
//             return this->solver(arg);
//         }
//     }
// }