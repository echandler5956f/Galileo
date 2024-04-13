#include <galileo/opt/PseudospectralSegment.h>

namespace galileo
{
    namespace opt
    {
        PseudospectralSegment::PseudospectralSegment(std::shared_ptr<GeneralProblemData> problem, casadi::Function F, casadi::Function L, std::shared_ptr<States> st_m, int d, int knot_num, double h)
        {
            auto Fint = problem->Fint;
            auto Fdiff = problem->Fdiff;

            assert(d > 0 && d < 10 && "d must be greater than 0 and less than 10");
            assert(h > 0 && "h must be a positive duration");

            assert(F.n_in() == 2 && "F must have 2 inputs");
            assert(F.n_out() == 1 && "F must have 1 output");
            assert(Fint.n_in() == 3 && "Fint must have 3 inputs");
            assert(Fint.n_out() == 1 && "Fint must have 1 output");
            assert(Fdiff.n_in() == 3 && "Fdiff must have 3 inputs");
            assert(Fdiff.n_out() == 1 && "Fdiff must have 1 output");
            assert(L.n_in() == 2 && "L must have 2 inputs");
            assert(L.n_out() == 1 && "L must have 1 output");

            F.assert_size_in(0, st_m->nx, 1);
            F.assert_size_in(1, st_m->nu, 1);
            F.assert_size_out(0, st_m->ndx, 1);

            L.assert_size_in(0, st_m->nx, 1);
            L.assert_size_in(1, st_m->nu, 1);
            L.assert_size_out(0, 1, 1);

            Fint.assert_size_in(0, st_m->nx, 1);
            Fint.assert_size_in(1, st_m->ndx, 1);
            Fint.assert_size_in(2, 1, 1);
            Fint.assert_size_out(0, st_m->nx, 1);

            Fdiff.assert_size_in(0, st_m->nx, 1);
            Fdiff.assert_size_in(1, st_m->nx, 1);
            Fdiff.assert_size_in(2, 1, 1);
            Fdiff.assert_size_out(0, st_m->ndx, 1);

            this->knot_num_ = knot_num;
            this->h_ = h;
            this->st_m_ = st_m;
            this->Fint_ = Fint;
            this->Fdiff_ = Fdiff;
            this->F_ = F;
            this->L_ = L;
            this->T_ = (knot_num)*h;

            InitializeExpressionVariables(d);
        }

        void PseudospectralSegment::InitializeExpressionVariables(int d)
        {
            expr_v_.dXc.clear();
            expr_v_.Uc.clear();

            dX_poly_ = LagrangePolynomial(d);
            U_poly_ = LagrangePolynomial(d - 1);

            for (int j = 0; j < dX_poly_.d; ++j)
            {
                expr_v_.dXc.push_back(casadi::SX::sym("dXc_" + std::to_string(j), st_m_->ndx, 1));
                if (j < U_poly_.d)
                {
                    expr_v_.Uc.push_back(casadi::SX::sym("Uc_" + std::to_string(j), st_m_->nu, 1));
                }
            }
            expr_v_.dX0 = casadi::SX::sym("dX0", st_m_->ndx, 1);
            expr_v_.X0 = casadi::SX::sym("X0", st_m_->nx, 1);
            expr_v_.U0 = casadi::SX::sym("U0", st_m_->nu, 1);
            expr_v_.Lc = casadi::SX::sym("Lc", 1, 1);
        }

        void PseudospectralSegment::InitializeSegmentTimeVector(casadi::DM &trajectory_times)
        {
            std::vector<double> vec_knot_times;
            std::vector<double> vec_collocation_times;
            std::vector<double> vec_segment_times;
            double kh = 0;
            for (int k = 0; k < knot_num_; ++k)
            {
                kh = k * h_;
                vec_knot_times.push_back(kh);
                for (int i = 0; i < dX_poly_.d; ++i)
                {
                    vec_collocation_times.push_back(kh + dX_poly_.tau_root[i + 1] * h_);
                }
            }
            vec_knot_times.push_back(T_);

            vec_segment_times.reserve(vec_knot_times.size() - 1 + vec_collocation_times.size());
            vec_segment_times.insert(vec_segment_times.end(), vec_knot_times.begin(), vec_knot_times.end() - 1);
            vec_segment_times.insert(vec_segment_times.end(), vec_collocation_times.begin(), vec_collocation_times.end());
            std::sort(vec_segment_times.begin(), vec_segment_times.end());

            tools::vectorToCasadi<casadi::DM>(vec_segment_times, (dX_poly_.d + 1) * knot_num_, 1, dXtimes_.segment_times);
            tools::vectorToCasadi<casadi::DM>(vec_collocation_times, dX_poly_.d * knot_num_, 1, dXtimes_.collocation_times);
            tools::vectorToCasadi<casadi::DM>(vec_knot_times, knot_num_ + 1, 1, dXtimes_.knot_times);

            double start_time = 0.0;
            if (trajectory_times.is_empty() == false)
            {
                start_time = trajectory_times(trajectory_times.size1() - 1, 0).get_elements()[0];
                dXtimes_.segment_times += start_time;
                trajectory_times = vertcat(trajectory_times, dXtimes_.segment_times);
            }
            else
            {
                trajectory_times = dXtimes_.segment_times;
            }
            dXtimes_.collocation_times += start_time;
            dXtimes_.knot_times += start_time;
        }

        void PseudospectralSegment::InitializeInputTimeVector(casadi::DM &trajectory_times)
        {
            std::vector<double> vec_u_knot_times;
            std::vector<double> vec_u_collocation_times;
            std::vector<double> vec_u_segment_times;
            double kh = 0;
            for (int k = 0; k < knot_num_; ++k)
            {
                kh = k * h_;
                vec_u_knot_times.push_back(kh);
                for (int i = 0; i < U_poly_.d; ++i)
                {
                    vec_u_collocation_times.push_back(kh + U_poly_.tau_root[i + 1] * h_);
                }
            }
            vec_u_knot_times.push_back(T_);

            vec_u_segment_times.reserve(vec_u_knot_times.size() - 1 + vec_u_collocation_times.size());
            vec_u_segment_times.insert(vec_u_segment_times.end(), vec_u_knot_times.begin(), vec_u_knot_times.end() - 1);
            vec_u_segment_times.insert(vec_u_segment_times.end(), vec_u_collocation_times.begin(), vec_u_collocation_times.end());
            std::sort(vec_u_segment_times.begin(), vec_u_segment_times.end());

            tools::vectorToCasadi<casadi::DM>(vec_u_segment_times, (U_poly_.d + 1) * knot_num_, 1, Utimes_.segment_times);
            tools::vectorToCasadi<casadi::DM>(vec_u_collocation_times, U_poly_.d * knot_num_, 1, Utimes_.collocation_times);
            tools::vectorToCasadi<casadi::DM>(vec_u_knot_times, knot_num_ + 1, 1, Utimes_.knot_times);

            double start_time = dXtimes_.knot_times(0, 0).get_elements()[0];
            Utimes_.segment_times += start_time;
            Utimes_.collocation_times += start_time;
            Utimes_.knot_times += start_time;
        }

        void PseudospectralSegment::InitializeKnotSegments(casadi::DM x0_global)
        {
            assert(x0_global.size1() == st_m_->nx && x0_global.size2() == 1 && "x0 must be a column std::vector of size nx");
            x0_global_ = x0_global;

            ps_vars_.dXc_vec.clear();
            ps_vars_.Uc_vec.clear();
            ps_vars_.dX0_vec.clear();
            ps_vars_.X0_vec.clear();
            ps_vars_.U0_vec.clear();
            for (int k = 0; k < knot_num_; ++k)
            {
                ps_vars_.dXc_vec.push_back(casadi::MX::sym("dXc_" + std::to_string(k), st_m_->ndx * dX_poly_.d, 1));
                ps_vars_.Uc_vec.push_back(casadi::MX::sym("U_" + std::to_string(k), st_m_->nu * U_poly_.d, 1));
            }

            /*We do knot_num + 1 so we have a decision variable for the final state. knot_num -1 is the number of knot segments, which corresponds to knot_num knot points*/
            for (int k = 0; k < knot_num_ + 1; ++k)
            {
                ps_vars_.dX0_vec.push_back(casadi::MX::sym("dX0_" + std::to_string(k), st_m_->ndx, 1));
                ps_vars_.X0_vec.push_back(Fint_(casadi::MXVector{x0_global_, ps_vars_.dX0_vec[k], 1.0}).at(0));
                ps_vars_.U0_vec.push_back(casadi::MX::sym("U0_" + std::to_string(k), st_m_->nu, 1));
            }
        }

        void PseudospectralSegment::InitializeExpressionGraph(std::vector<ConstraintData> G, std::shared_ptr<DecisionData> Wdata)
        {
            /*Collocation equations*/
            casadi::SXVector eq;
            /*State at the end of the collocation interval*/
            casadi::SX dXf = dX_poly_.D[0] * expr_v_.dX0;
            casadi::SX uf = U_poly_.D[0] * expr_v_.U0;
            /*Cost at the end of the collocation interval*/
            casadi::SX Qf = 0;
            /*Actual state at collocation points*/
            casadi::SXVector x_at_c;
            /*U interpolated at the dx polynomial collocation points*/
            casadi::SXVector u_at_c;
            casadi::SXVector tmp_x;
            casadi::SXVector tmp_dx;
            casadi::SXVector tmp_u = {expr_v_.U0};
            tmp_x.push_back(expr_v_.X0);
            tmp_dx.push_back(expr_v_.dX0);
            for (int j = 0; j < U_poly_.d; ++j)
            {
                tmp_u.push_back(expr_v_.Uc[j]);
                uf += U_poly_.D[j + 1] * expr_v_.Uc[j];
            }

            for (int j = 0; j < dX_poly_.d; ++j)
            {
                double dt_j = (dX_poly_.tau_root[j + 1] - dX_poly_.tau_root[j]) * h_;
                /*Expression for the state derivative at the collocation point*/
                casadi::SX dxp = dX_poly_.C[0][j + 1] * expr_v_.dX0;
                for (int r = 0; r < dX_poly_.d; ++r)
                {
                    dxp += dX_poly_.C[r + 1][j + 1] * expr_v_.dXc[r];
                }
                /*dXc must exist in a Euclidean space, but we need x_c in order to evaluate the objective. Fint can simply return dXc[j] if the states are already Euclidean*/
                casadi::SX x_c = Fint_(casadi::SXVector{expr_v_.X0, expr_v_.dXc[j], dt_j}).at(0);
                casadi::SX u_c = U_poly_.BarycentricInterpolation(dX_poly_.tau_root[j], tmp_u);

                x_at_c.push_back(x_c);
                u_at_c.push_back(u_c);
                tmp_x.push_back(x_c);
                tmp_dx.push_back(expr_v_.dXc[j]);

                /*Append collocation equations*/
                eq.push_back(h_ * F_(casadi::SXVector{x_c, u_c}).at(0) - dxp);

                /*Add cost contribution*/
                casadi::SXVector L_out = L_(casadi::SXVector{x_c, u_c});
                /*This is fine as long as the cost is not related to the Lie Group elements. See the state integrator and dX for clarity*/
                Qf += dX_poly_.B[j + 1] * L_out.at(0) * h_;

                dXf += dX_poly_.D[j + 1] * expr_v_.dXc[j];
            }

            casadi::SX vcat_dXc = vertcat(expr_v_.dXc);
            casadi::SX vcat_Uc = vertcat(expr_v_.Uc);
            casadi::SXVector function_inputs = {expr_v_.X0, vcat_dXc, expr_v_.dX0, expr_v_.U0, vcat_Uc};

            casadi::Function collocation_constraint = casadi::Function("feq",
                                                                       function_inputs,
                                                                       casadi::SXVector{vertcat(eq)});

            casadi::Function xf_constraint = casadi::Function("fxf",
                                                              function_inputs,
                                                              casadi::SXVector{dXf});

            casadi::Function uf_constraint = casadi::Function("fuf",
                                                              function_inputs,
                                                              casadi::SXVector{uf});

            casadi::Function q_cost = casadi::Function("fxq", casadi::SXVector{expr_v_.Lc, expr_v_.X0, vcat_dXc, expr_v_.dX0, expr_v_.U0, vcat_Uc},
                                                       casadi::SXVector{expr_v_.Lc + Qf});

            /*Implicit discrete-time equations*/
            ps_funcs_.col_con_map = collocation_constraint.map(knot_num_, "openmp");
            /*When you evaluate this map, subtract by the knot points list offset by 1 to be correct*/
            ps_funcs_.xf_con_map = xf_constraint.map(knot_num_, "openmp");
            ps_funcs_.uf_con_map = uf_constraint.map(knot_num_, "openmp");
            ps_funcs_.q_cost_fold = q_cost.fold(knot_num_);

            sol_map_func_ = casadi::Function("sol_map",
                                             function_inputs,
                                             casadi::SXVector{horzcat(tmp_x), horzcat(tmp_u)})
                                .map(knot_num_, "serial");

            casadi_int N = ps_funcs_.col_con_map.size1_out(0) * ps_funcs_.col_con_map.size2_out(0) +
                           ps_funcs_.xf_con_map.size1_out(0) * ps_funcs_.xf_con_map.size2_out(0) +
                           ps_funcs_.uf_con_map.size1_out(0) * ps_funcs_.uf_con_map.size2_out(0);
            casadi_int tmp = N;

            std::vector<tuple_size_t> ranges_G;

            casadi::SXVector tmap_symbolic_input = casadi::SXVector{horzcat(x_at_c), horzcat(u_at_c)};
            /*Map the constraint to each collocation point, and then map the mapped constraint to each knot segment*/
            for (size_t i = 0; i < G.size(); ++i)
            {
                ConstraintData g_data = G[i];

                assert(g_data.G.n_in() == 2 && "G must have 2 inputs");
                g_data.G.assert_size_in(0, st_m_->nx, 1);
                g_data.G.assert_size_in(1, st_m_->nu, 1);
                assert(g_data.lower_bound.n_in() == 1 && "G lower_bound must have 1 inputs");
                assert(g_data.lower_bound.n_out() == 1 && "G lower_bound must have 1 output");
                g_data.lower_bound.assert_size_in(0, 1, 1);
                casadi::Function tmap = casadi::Function(g_data.G.name() + "_map",
                                                         function_inputs,
                                                         casadi::SXVector{vertcat(g_data.G.map(dX_poly_.d, "serial")((tmap_symbolic_input)))})
                                            .map(knot_num_, "serial");
                ps_funcs_.gcon_maps.push_back(tmap);
                ranges_G.push_back(tuple_size_t(N, N + tmap.size1_out(0) * tmap.size2_out(0)));
                N += tmap.size1_out(0) * tmap.size2_out(0);
            }

            nlp_data_.lbg.push_back(casadi::DM::zeros(tmp, 1));
            nlp_data_.ubg.push_back(casadi::DM::zeros(tmp, 1));

            ConstraintData g_data;

            for (std::size_t i = 0; i < G.size(); ++i)
            {
                g_data = G[i];
                nlp_data_.lbg.push_back(casadi::DM::reshape(vertcat(g_data.lower_bound.map(knot_num_ * (dX_poly_.d), "serial")(dXtimes_.collocation_times)), std::get<1>(ranges_G[i]) - std::get<0>(ranges_G[i]), 1));
                nlp_data_.ubg.push_back(casadi::DM::reshape(vertcat(g_data.upper_bound.map(knot_num_ * (dX_poly_.d), "serial")(dXtimes_.collocation_times)), std::get<1>(ranges_G[i]) - std::get<0>(ranges_G[i]), 1));
            }

            int Ndxknot = st_m_->ndx * (knot_num_ + 1);
            int Ndx = st_m_->ndx * (dX_poly_.d + 1) * knot_num_ + st_m_->ndx;
            int Ndxcol = Ndx - Ndxknot;

            int Nuknot = st_m_->nu * (knot_num_ + 1);
            int Nu = st_m_->nu * (U_poly_.d + 1) * knot_num_ + st_m_->nu;
            int Nucol = Nu - Nuknot;

            /*Transform initial guess for x to an initial guess for dx, using f_diff, the inverse of f_int*/
            casadi::MX xkg_sym = casadi::MX::sym("xkg", st_m_->nx, 1);
            casadi::MX xckg_sym = casadi::MX::sym("Xckg", st_m_->nx * dX_poly_.d, 1);
            if (!Wdata->initial_guess.is_null())
            {
                casadi::DM xg = Wdata->initial_guess.map(knot_num_ + 1, "serial")(dXtimes_.knot_times).at(0);
                casadi::Function dxg_func = casadi::Function("xg_fun", casadi::MXVector{xkg_sym}, casadi::MXVector{Fdiff_(casadi::MXVector{x0_global_, xkg_sym, 1.0}).at(0)})
                                                .map(knot_num_ + 1, "serial");
                nlp_data_.w0.push_back(casadi::DM::reshape(dxg_func(casadi::DMVector{xg}).at(0), Ndxknot, 1));
                /*The transformation of xc to dxc is a slightly less trivial. While x_k = fint(x0_init, dx_k), for xc_k, we have xc_k = fint(x_k, dxc_k) which is equivalent to xc_k = fint(fint(x0_init, dx_k), dxc_k).
                Thus, dxc_k = fdiff(fint(x0_init, dx_k), xc_k)). This could be done with maps like above, but it is not necessary.*/
                casadi::DM xc_g = Wdata->initial_guess.map((dX_poly_.d) * knot_num_, "serial")(dXtimes_.collocation_times).at(0);
                for (casadi_int i = 0; i < knot_num_; ++i)
                {
                    casadi::DM xk = xg(casadi::Slice(i * st_m_->nx, (i + 1) * st_m_->nx));
                    casadi::DM xck = xc_g(casadi::Slice(i * st_m_->nx * dX_poly_.d, (i + 1) * st_m_->nx * dX_poly_.d));
                    for (casadi_int j = 0; j < dX_poly_.d; ++j)
                    {
                        nlp_data_.w0.push_back(reshape(Fdiff_(casadi::DMVector{xk, xck(casadi::Slice(j * st_m_->nx, (j + 1) * st_m_->nx)), h_}).at(0), st_m_->ndx, 1));
                    }
                }

                nlp_data_.w0.push_back(casadi::DM::reshape(Wdata->initial_guess.map(knot_num_ + 1, "serial")(Utimes_.knot_times).at(1), Nuknot, 1));
                nlp_data_.w0.push_back(casadi::DM::reshape(Wdata->initial_guess.map((U_poly_.d) * knot_num_, "serial")(Utimes_.collocation_times).at(1), Nucol, 1));
            }
            else
            {
                nlp_data_.w0.push_back(casadi::DM::zeros(Ndx + Nu, 1));
            }

            if (!Wdata->lower_bound.is_null() && !Wdata->upper_bound.is_null())
            {
                nlp_data_.lbw.push_back(casadi::DM::reshape(Wdata->lower_bound.map(knot_num_ + 1, "serial")(dXtimes_.knot_times).at(0), Ndxknot, 1));
                nlp_data_.ubw.push_back(casadi::DM::reshape(Wdata->upper_bound.map(knot_num_ + 1, "serial")(dXtimes_.knot_times).at(0), Ndxknot, 1));
                nlp_data_.lbw.push_back(casadi::DM::reshape(Wdata->lower_bound.map((dX_poly_.d) * knot_num_, "serial")(dXtimes_.collocation_times).at(0), Ndxcol, 1));
                nlp_data_.ubw.push_back(casadi::DM::reshape(Wdata->upper_bound.map((dX_poly_.d) * knot_num_, "serial")(dXtimes_.collocation_times).at(0), Ndxcol, 1));

                nlp_data_.lbw.push_back(casadi::DM::reshape(Wdata->lower_bound.map(knot_num_ + 1, "serial")(Utimes_.knot_times).at(1), Nuknot, 1));
                nlp_data_.ubw.push_back(casadi::DM::reshape(Wdata->upper_bound.map(knot_num_ + 1, "serial")(Utimes_.knot_times).at(1), Nuknot, 1));
                nlp_data_.lbw.push_back(casadi::DM::reshape(Wdata->lower_bound.map((U_poly_.d) * knot_num_, "serial")(Utimes_.collocation_times).at(1), Nucol, 1));
                nlp_data_.ubw.push_back(casadi::DM::reshape(Wdata->upper_bound.map((U_poly_.d) * knot_num_, "serial")(Utimes_.collocation_times).at(1), Nucol, 1));
            }
            else
            {
                nlp_data_.lbw.push_back(-casadi::DM::inf(Ndx + Nu, 1));
                nlp_data_.ubw.push_back(casadi::DM::inf(Ndx + Nu, 1));
            }
        }

        void PseudospectralSegment::EvaluateExpressionGraph()
        {
            casadi::MXVector result;

            casadi::MX xs = ProcessVector(ps_vars_.X0_vec);
            casadi::MX us = ProcessVector(ps_vars_.U0_vec);
            casadi::MX dxs = ProcessVector(ps_vars_.dX0_vec);
            casadi::MX dxcs = horzcat(ps_vars_.dXc_vec);
            casadi::MX ucs = horzcat(ps_vars_.Uc_vec);
            casadi::MX xs_offset = ProcessOffsetVector(ps_vars_.X0_vec);
            casadi::MX dxs_offset = ProcessOffsetVector(ps_vars_.dX0_vec);
            casadi::MX us_offset = ProcessOffsetVector(ps_vars_.U0_vec);

            casadi::MXVector solmap_result = sol_map_func_(casadi::MXVector{xs, dxcs, dxs, us, ucs});
            casadi::MX all_xs = solmap_result.at(0);
            casadi::MX all_us = solmap_result.at(1);

            /*This section cannot get much faster, it is bounded by the time to evaluate the constraint*/
            casadi::MX col_con_mat = ps_funcs_.col_con_map(casadi::MXVector{xs, dxcs, dxs, us, ucs}).at(0);
            casadi::MX xf_con_mat = ps_funcs_.xf_con_map(casadi::MXVector{xs, dxcs, dxs, us, ucs}).at(0);
            casadi::MX uf_con_mat = ps_funcs_.uf_con_map(casadi::MXVector{xs, dxcs, dxs, us, ucs}).at(0);
            dxs_offset = reshape(dxs_offset, dxs_offset.size1() * dxs_offset.size2(), 1);
            us_offset = reshape(us_offset, us_offset.size1() * us_offset.size2(), 1);

            result.push_back(reshape(col_con_mat, col_con_mat.size1() * col_con_mat.size2(), 1));
            result.push_back(reshape(xf_con_mat, xf_con_mat.size1() * xf_con_mat.size2(), 1) -
                             dxs_offset);
            result.push_back(reshape(uf_con_mat, uf_con_mat.size1() * uf_con_mat.size2(), 1) -
                             us_offset);

            for (size_t i = 0; i < ps_funcs_.gcon_maps.size(); ++i)
            {
                casadi::MX g_con_mat = ps_funcs_.gcon_maps[i](casadi::MXVector{xs, dxcs, dxs, us, ucs}).at(0);
                result.push_back(reshape(g_con_mat, g_con_mat.size1() * g_con_mat.size2(), 1));
            }

            nlp_data_.J = ps_funcs_.q_cost_fold(casadi::MXVector{0., xs, dxcs, dxs, us, ucs}).at(0);

            nlp_data_.g.insert(nlp_data_.g.end(), result.begin(), result.end());

            nlp_data_.w.insert(nlp_data_.w.end(), ps_vars_.dX0_vec.begin(), ps_vars_.dX0_vec.end());
            nlp_data_.w.insert(nlp_data_.w.end(), ps_vars_.dXc_vec.begin(), ps_vars_.dXc_vec.end());
            nlp_data_.w.insert(nlp_data_.w.end(), ps_vars_.U0_vec.begin(), ps_vars_.U0_vec.end());
            nlp_data_.w.insert(nlp_data_.w.end(), ps_vars_.Uc_vec.begin(), ps_vars_.Uc_vec.end());

            std::cout << "Pseudospectral w size: " << vertcat(nlp_data_.w).size1() << std::endl;
            get_sol_func_ = casadi::Function("func",
                                             casadi::MXVector({vertcat(nlp_data_.w)}),
                                             casadi::MXVector({all_xs, all_us}));
        }

        casadi::MXVector PseudospectralSegment::ExtractSolution(casadi::MX &w) const
        {
            return get_sol_func_(casadi::MXVector{w});
        }

        void PseudospectralSegment::FillNLPData(NLPData &nlp_data)
        {
            /*Where w of this segment starts*/
            size_t w_it_size = nlp_data.w.size();
            size_t w_size = 0;
            /*Accumulate the size of each symbolic element in w*/
            if (w_it_size != 0)
                w_size = size_t(accumulate(nlp_data.w.begin(), nlp_data.w.end(), 0.0, [](int sum, const casadi::MX &item)
                                           { return sum + item.size1() * item.size2(); }));

            nlp_data.w.insert(nlp_data.w.end(), nlp_data_.w.begin(), nlp_data_.w.end());
            w_range_ = tuple_size_t(w_size, accumulate(nlp_data.w.begin(), nlp_data.w.end(), 0.0, [](int sum, const casadi::MX &item)
                                                       { return sum + item.size1() * item.size2(); }));

            nlp_data.lbw.insert(nlp_data.lbw.end(), make_move_iterator(nlp_data_.lbw.begin()), make_move_iterator(nlp_data_.lbw.end()));
            nlp_data.ubw.insert(nlp_data.ubw.end(), make_move_iterator(nlp_data_.ubw.begin()), make_move_iterator(nlp_data_.ubw.end()));
            nlp_data.w0.insert(nlp_data.w0.end(), make_move_iterator(nlp_data_.w0.begin()), make_move_iterator(nlp_data_.w0.end()));
            nlp_data.g.insert(nlp_data.g.end(), make_move_iterator(nlp_data_.g.begin()), make_move_iterator(nlp_data_.g.end()));
            nlp_data.lbg.insert(nlp_data.lbg.end(), make_move_iterator(nlp_data_.lbg.begin()), make_move_iterator(nlp_data_.lbg.end()));
            nlp_data.ubg.insert(nlp_data.ubg.end(), make_move_iterator(nlp_data_.ubg.begin()), make_move_iterator(nlp_data_.ubg.end()));
            nlp_data.J += nlp_data_.J;
        }

        casadi::MX PseudospectralSegment::ProcessVector(casadi::MXVector &vec) const
        {
            casadi::MXVector temp = vec;
            temp.pop_back();
            return horzcat(temp);
        }

        casadi::MX PseudospectralSegment::ProcessOffsetVector(casadi::MXVector &vec) const
        {
            casadi::MXVector temp = vec;
            temp.erase(temp.begin());
            return horzcat(temp);
        }
    }
}