#include <galileo/opt/PseudospectralSegment.h>

namespace galileo
{
    namespace opt
    {
        PseudospectralSegment::PseudospectralSegment(std::shared_ptr<GeneralProblemData> problem, casadi::Function F_, std::shared_ptr<States> st_m_, int d, int knot_num_, double h_)
        {
            auto Fint_ = problem->Fint;
            auto Fdif_ = problem->Fdif;
            auto L_ = problem->L;

            assert(d > 0 && d < 10 && "d must be greater than 0 and less than 10");
            assert(h_ > 0 && "h must be a positive duration");

            assert(F_.n_in() == 2 && "F must have 2 inputs");
            assert(F_.n_out() == 1 && "F must have 1 output");
            assert(Fint_.n_in() == 3 && "Fint must have 3 inputs");
            assert(Fint_.n_out() == 1 && "Fint must have 1 output");
            assert(Fdif_.n_in() == 3 && "Fdif must have 3 inputs");
            assert(Fdif_.n_out() == 1 && "Fdif must have 1 output");
            assert(L_.n_in() == 2 && "L must have 2 inputs");
            assert(L_.n_out() == 1 && "L must have 1 output");

            F_.assert_size_in(0, st_m_->nx, 1);
            F_.assert_size_in(1, st_m_->nu, 1);
            F_.assert_size_out(0, st_m_->ndx, 1);

            L_.assert_size_in(0, st_m_->nx, 1);
            L_.assert_size_in(1, st_m_->nu, 1);
            L_.assert_size_out(0, 1, 1);

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

            initializeExpressionVariables(d);
        }

        void PseudospectralSegment::initializeExpressionVariables(int d)
        {
            dXc.clear();
            Uc.clear();

            dX_poly = LagrangePolynomial(d);
            U_poly = LagrangePolynomial(0);

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

        void PseudospectralSegment::initializeLocalTimeVector(casadi::DM &global_times)
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
            initializeInputTimeVector();
        }

        void PseudospectralSegment::initializeInputTimeVector()
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

        void PseudospectralSegment::initializeKnotSegments(casadi::DM x0_global_, casadi::MX x0_local_)
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
                X0_var_vec.push_back(Fint(casadi::MXVector{x0_global_, dX0_var_vec[k], 1.0}).at(0));
            }
        }

        void PseudospectralSegment::initializeExpressionGraph(std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu)
        {
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
                assert(g_data->lower_bound.n_in() == 1 && "G lower_bound must have 1 inputs");
                assert(g_data->lower_bound.n_out() == 1 && "G lower_bound must have 1 output");
                g_data->lower_bound.assert_size_in(0, 1, 1);
                auto tmap = casadi::Function(g_data->G.name() + "_map",
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

            for (std::size_t i = 0; i < G.size(); ++i)
            {
                auto g_data = G[i];
                general_lbg(casadi::Slice(casadi_int(std::get<0>(ranges_G[i])), casadi_int(std::get<1>(ranges_G[i]))), 0) =
                    casadi::DM::reshape(vertcat(g_data->lower_bound.map(knot_num * (dX_poly.d), "serial")(collocation_times)), std::get<1>(ranges_G[i]) - std::get<0>(ranges_G[i]), 1);
                general_ubg(casadi::Slice(casadi_int(std::get<0>(ranges_G[i])), casadi_int(std::get<1>(ranges_G[i]))), 0) =
                    casadi::DM::reshape(vertcat(g_data->upper_bound.map(knot_num * (dX_poly.d), "serial")(collocation_times)), std::get<1>(ranges_G[i]) - std::get<0>(ranges_G[i]), 1);
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
                w0(casadi::Slice(0, Ndxknot)) = casadi::DM::reshape(dxg_func(casadi::DMVector{xg}).at(0), Ndxknot, 1);
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
                general_lbw(casadi::Slice(0, Ndxknot)) = casadi::DM::reshape(vertcat(Wx->lower_bound.map(knot_num + 1, "serial")(knot_times)), Ndxknot, 1);
                general_ubw(casadi::Slice(0, Ndxknot)) = casadi::DM::reshape(vertcat(Wx->upper_bound.map(knot_num + 1, "serial")(knot_times)), Ndxknot, 1);
                general_lbw(casadi::Slice(Ndxknot, Ndx)) = casadi::DM::reshape(vertcat(Wx->lower_bound.map((dX_poly.d) * knot_num, "serial")(collocation_times)), Ndxcol, 1);
                general_ubw(casadi::Slice(Ndxknot, Ndx)) = casadi::DM::reshape(vertcat(Wx->upper_bound.map((dX_poly.d) * knot_num, "serial")(collocation_times)), Ndxcol, 1);
            }

            if (!Wu->initial_guess.is_null())
            {
                w0(casadi::Slice(Ndx, Ndx + Nu)) = casadi::DM::reshape(vertcat(Wu->initial_guess.map(U_poly.d * knot_num, "serial")(u_times)), Nu, 1);
            }

            if (!Wu->lower_bound.is_null() && !Wu->upper_bound.is_null())
            {
                general_lbw(casadi::Slice(Ndx, Ndx + Nu)) = casadi::DM::reshape(vertcat(Wu->lower_bound.map(U_poly.d * knot_num, "serial")(u_times)), Nu, 1);
                general_ubw(casadi::Slice(Ndx, Ndx + Nu)) = casadi::DM::reshape(vertcat(Wu->upper_bound.map(U_poly.d * knot_num, "serial")(u_times)), Nu, 1);
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

        void PseudospectralSegment::evaluateExpressionGraph(casadi::MX &J0, casadi::MXVector &w, casadi::MXVector &g)
        {
            assert(J0.size1() == 1 && J0.size2() == 1 && "J0 must be a scalar");

            casadi::MXVector result;

            casadi::MX xs = processVector(X0_var_vec);
            casadi::MX dxs = processVector(dX0_var_vec);
            casadi::MX dxcs = horzcat(dXc_var_vec);
            casadi::MX us = horzcat(U_var_vec);
            casadi::MX xs_offset = processOffsetVector(X0_var_vec);
            casadi::MX dxs_offset = processOffsetVector(dX0_var_vec);

            casadi::MXVector solmap_result = sol_map_func(casadi::MXVector{xs, dxcs, dxs, us});
            casadi::MX all_xs = solmap_result.at(0);
            casadi::MX all_us = solmap_result.at(1);

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
            size_t w_it_size = w.size();
            size_t w_size = 0;
            /*Accumulate the size of each symbolic element in w*/
            if (w_it_size != 0)
                w_size = size_t(accumulate(w.begin(), w.end(), 0.0, [](int sum, const casadi::MX &item)
                                    { return sum + item.size1() * item.size2(); }));

            /*Use move to avoid copying the vectors. Reserve space for w in advance outside of PseudospectralSegment.*/
            w.insert(w.end(), make_move_iterator(dX0_var_vec.begin()), make_move_iterator(dX0_var_vec.end()));
            w.insert(w.end(), make_move_iterator(dXc_var_vec.begin()), make_move_iterator(dXc_var_vec.end()));
            w.insert(w.end(), make_move_iterator(U_var_vec.begin()), make_move_iterator(U_var_vec.end()));
                                                  
            w_range = tuple_size_t(w_size, accumulate(w.begin(), w.end(), 0.0, [](int sum, const casadi::MX &item)
                                                      { return sum + item.size1() * item.size2(); }));
            get_sol_func = casadi::Function("func",
                                            casadi::MXVector({vertcat(w)}),
                                            casadi::MXVector({all_xs, all_us}));
        }

        casadi::MXVector PseudospectralSegment::extractSolution(casadi::MX &w) const
        {
            return get_sol_func(casadi::MXVector{w});
        }

        casadi::MX PseudospectralSegment::getInitialStateDeviant() const
        {
            return dX0_var_vec.front();
        }

        casadi::MX PseudospectralSegment::getInitialState() const
        {
            return X0_var_vec.front();
        }

        casadi::MX PseudospectralSegment::getFinalStateDeviant() const
        {
            return dX0_var_vec.back();
        }

        casadi::MX PseudospectralSegment::getFinalState() const
        {
            return X0_var_vec.back();
        }

        casadi::DM PseudospectralSegment::getLocalTimes() const
        {
            return local_times;
        }

        casadi::DM PseudospectralSegment::getInputTimes() const
        {
            return u_times;
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