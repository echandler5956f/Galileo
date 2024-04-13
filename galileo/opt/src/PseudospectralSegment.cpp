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

            initializeExpressionVariables(d);
        }

        void PseudospectralSegment::initializeExpressionVariables(int d)
        {
            dXc_.clear();
            Uc_.clear();

            dX_poly_ = LagrangePolynomial(d);
            U_poly_ = LagrangePolynomial(d - 1);

            for (int j = 0; j < dX_poly_.d; ++j)
            {
                dXc_.push_back(casadi::SX::sym("dXc_" + std::to_string(j), st_m_->ndx, 1));
                if (j < U_poly_.d)
                {
                    Uc_.push_back(casadi::SX::sym("Uc_" + std::to_string(j), st_m_->nu, 1));
                }
            }
            dX0_ = casadi::SX::sym("dX0", st_m_->ndx, 1);
            X0_ = casadi::SX::sym("X0", st_m_->nx, 1);
            U0_ = casadi::SX::sym("U0", st_m_->nu, 1);
            Lc_ = casadi::SX::sym("Lc", 1, 1);
        }

        void PseudospectralSegment::initializeSegmentTimeVector(casadi::DM &trajectory_times)
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

            tools::vectorToCasadi<casadi::DM>(vec_segment_times, (dX_poly_.d + 1) * knot_num_, 1, segment_times_);
            tools::vectorToCasadi<casadi::DM>(vec_collocation_times, dX_poly_.d * knot_num_, 1, collocation_times_);
            tools::vectorToCasadi<casadi::DM>(vec_knot_times, knot_num_ + 1, 1, knot_times_);

            double start_time = 0.0;
            if (trajectory_times.is_empty() == false)
            {
                start_time = trajectory_times(trajectory_times.size1() - 1, 0).get_elements()[0];
                segment_times_ += start_time;
                trajectory_times = vertcat(trajectory_times, segment_times_);
            }
            else
            {
                trajectory_times = segment_times_;
            }
            t_range_ = tuple_size_t(trajectory_times.size1() - segment_times_.size1(), trajectory_times.size1());
            collocation_times_ += start_time;
            knot_times_ += start_time;
        }

        // This has a high chance of being different than initializeSegmentTimeVector in the future
        void PseudospectralSegment::initializeInputTimeVector(casadi::DM &trajectory_times)
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

            tools::vectorToCasadi<casadi::DM>(vec_u_segment_times, (U_poly_.d + 1) * knot_num_, 1, u_segment_times_);
            tools::vectorToCasadi<casadi::DM>(vec_u_collocation_times, U_poly_.d * knot_num_, 1, u_collocation_times_);
            tools::vectorToCasadi<casadi::DM>(vec_u_knot_times, knot_num_ + 1, 1, u_knot_times_);

            double start_time = knot_times_(0, 0).get_elements()[0];
            u_segment_times_ += start_time;
            u_collocation_times_ += start_time;
            u_knot_times_ += start_time;
        }

        void PseudospectralSegment::initializeKnotSegments(casadi::DM x0_global)
        {
            assert(x0_global.size1() == st_m_->nx && x0_global.size2() == 1 && "x0 must be a column std::vector of size nx");
            x0_global_ = x0_global;

            dXc_var_vec_.clear();
            Uc_var_vec_.clear();
            dX0_var_vec_.clear();
            X0_var_vec_.clear();
            U0_var_vec_.clear();
            U_at_c_vec_.clear();
            x_at_c_vec_.clear();
            for (int k = 0; k < knot_num_; ++k)
            {
                dXc_var_vec_.push_back(casadi::MX::sym("dXc_" + std::to_string(k), st_m_->ndx * dX_poly_.d, 1));
                Uc_var_vec_.push_back(casadi::MX::sym("U_" + std::to_string(k), st_m_->nu * U_poly_.d, 1));
            }

            /*We do knot_num + 1 so we have a decision variable for the final state. knot_num -1 is the number of knot segments, which corresponds to knot_num knot points*/
            for (int k = 0; k < knot_num_ + 1; ++k)
            {
                dX0_var_vec_.push_back(casadi::MX::sym("dX0_" + std::to_string(k), st_m_->ndx, 1));
                X0_var_vec_.push_back(Fint_(casadi::MXVector{x0_global_, dX0_var_vec_[k], 1.0}).at(0));
                U0_var_vec_.push_back(casadi::MX::sym("U0_" + std::to_string(k), st_m_->nu, 1));
            }
        }

        void PseudospectralSegment::initializeExpressionGraph(std::vector<ConstraintData> G, std::shared_ptr<DecisionData> Wdata)
        {
            /*Collocation equations*/
            casadi::SXVector eq;
            /*State at the end of the collocation interval*/
            casadi::SX dXf = dX_poly_.D[0] * dX0_;
            casadi::SX uf = U_poly_.D[0] * U0_;
            /*Cost at the end of the collocation interval*/
            casadi::SX Qf = 0;
            /*Actual state at collocation points*/
            casadi::SXVector x_at_c;
            /*U interpolated at the dx polynomial collocation points*/
            casadi::SXVector u_at_c;
            casadi::SXVector tmp_x;
            casadi::SXVector tmp_dx;
            casadi::SXVector tmp_u = {U0_};
            tmp_x.push_back(X0_);
            tmp_dx.push_back(dX0_);
            for (int j = 0; j < U_poly_.d; ++j)
            {
                tmp_u.push_back(Uc_[j]);
                uf += U_poly_.D[j + 1] * Uc_[j];
            }

            for (int j = 0; j < dX_poly_.d; ++j)
            {
                double dt_j = (dX_poly_.tau_root[j + 1] - dX_poly_.tau_root[j]) * h_;
                /*Expression for the state derivative at the collocation point*/
                casadi::SX dxp = dX_poly_.C[0][j + 1] * dX0_;
                for (int r = 0; r < dX_poly_.d; ++r)
                {
                    dxp += dX_poly_.C[r + 1][j + 1] * dXc_[r];
                }
                /*dXc must exist in a Euclidean space, but we need x_c in order to evaluate the objective. Fint can simply return dXc[j] if the states are already Euclidean*/
                casadi::SX x_c = Fint_(casadi::SXVector{X0_, dXc_[j], dt_j}).at(0);
                casadi::SX u_c = U_poly_.barycentricInterpolation(dX_poly_.tau_root[j], tmp_u);

                x_at_c.push_back(x_c);
                u_at_c.push_back(u_c);
                tmp_x.push_back(x_c);
                tmp_dx.push_back(dXc_[j]);

                /*Append collocation equations*/
                eq.push_back(h_ * F_(casadi::SXVector{x_c, u_c}).at(0) - dxp);

                /*Add cost contribution*/
                casadi::SXVector L_out = L_(casadi::SXVector{x_c, u_c});
                /*This is fine as long as the cost is not related to the Lie Group elements. See the state integrator and dX for clarity*/
                Qf += dX_poly_.B[j + 1] * L_out.at(0) * h_;

                dXf += dX_poly_.D[j + 1] * dXc_[j];
            }

            casadi::Dict opts;
            // // opts["cse"] = true;
            // opts["jit"] = true;
            // opts["jit_options.flags"] = "-Ofast -march=native -ffast-math";
            // opts["jit_options.compiler"] = "gcc";
            // // opts["jit_options.temp_suffix"] = false;
            // opts["compiler"] = "shell";

            casadi::SX vcat_dXc = vertcat(dXc_);
            casadi::SX vcat_Uc = vertcat(Uc_);
            casadi::SXVector function_inputs = {X0_, vcat_dXc, dX0_, U0_, vcat_Uc};


            casadi::Function collocation_constraint = casadi::Function("feq",
                                                           function_inputs,
                                                           casadi::SXVector{vertcat(eq)}, opts);

            casadi::Function xf_constraint = casadi::Function("fxf",
                                                  function_inputs,
                                                  casadi::SXVector{dXf}, opts);

            casadi::Function uf_constraint = casadi::Function("fuf",
                                                  function_inputs,
                                                  casadi::SXVector{uf}, opts);

            casadi::Function q_cost = casadi::Function("fxq", casadi::SXVector{Lc_, X0_, vcat_dXc, dX0_, U0_, vcat_Uc},
                                           casadi::SXVector{Lc_ + Qf}, opts);

            /*Implicit discrete-time equations*/
            collocation_constraint_map_ = collocation_constraint.map(knot_num_, "openmp");
            /*When you evaluate this map, subtract by the knot points list offset by 1 to be correct*/
            xf_constraint_map_ = xf_constraint.map(knot_num_, "openmp");
            uf_constraint_map_ = uf_constraint.map(knot_num_, "openmp");
            q_cost_fold_ = q_cost.fold(knot_num_);

            sol_map_func_ = casadi::Function("sol_map",
                                            function_inputs,
                                            casadi::SXVector{horzcat(tmp_x), horzcat(tmp_u)})
                               .map(knot_num_, "serial");

            casadi_int N = collocation_constraint_map_.size1_out(0) * collocation_constraint_map_.size2_out(0) +
                           xf_constraint_map_.size1_out(0) * xf_constraint_map_.size2_out(0) +
                           uf_constraint_map_.size1_out(0) * uf_constraint_map_.size2_out(0);
            casadi_int tmp = N;

            std::vector<tuple_size_t> ranges_G;

            casadi::SXVector tmap_symbolic_input = casadi::SXVector{horzcat(x_at_c), horzcat(u_at_c)};
            /*Map the constraint to each collocation point, and then map the mapped constraint to each knot segment*/
            for (size_t i = 0; i < G.size(); ++i)
            {
                ConstraintData  g_data = G[i];

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
                general_constraint_maps_.push_back(tmap);
                ranges_G.push_back(tuple_size_t(N, N + tmap.size1_out(0) * tmap.size2_out(0)));
                N += tmap.size1_out(0) * tmap.size2_out(0);
            }

            general_lbg_.resize(N, 1);
            general_ubg_.resize(N, 1);
            general_lbg_(casadi::Slice(0, tmp)) = casadi::DM::zeros(tmp, 1);
            general_ubg_(casadi::Slice(0, tmp)) = casadi::DM::zeros(tmp, 1);

            ConstraintData g_data;

            for (std::size_t i = 0; i < G.size(); ++i)
            {
                g_data = G[i];
                general_lbg_(casadi::Slice(casadi_int(std::get<0>(ranges_G[i])), casadi_int(std::get<1>(ranges_G[i]))), 0) =
                    casadi::DM::reshape(vertcat(g_data.lower_bound.map(knot_num_ * (dX_poly_.d), "serial")(collocation_times_)), std::get<1>(ranges_G[i]) - std::get<0>(ranges_G[i]), 1);
                general_ubg_(casadi::Slice(casadi_int(std::get<0>(ranges_G[i])), casadi_int(std::get<1>(ranges_G[i]))), 0) =
                    casadi::DM::reshape(vertcat(g_data.upper_bound.map(knot_num_ * (dX_poly_.d), "serial")(collocation_times_)), std::get<1>(ranges_G[i]) - std::get<0>(ranges_G[i]), 1);
            }

            int Ndxknot = st_m_->ndx * (knot_num_ + 1);
            int Ndx = st_m_->ndx * (dX_poly_.d + 1) * knot_num_ + st_m_->ndx;
            int Ndxcol = Ndx - Ndxknot;

            int Nuknot = st_m_->nu * (knot_num_ + 1);
            int Nu = st_m_->nu * (U_poly_.d + 1) * knot_num_ + st_m_->nu;
            int Nucol = Nu - Nuknot;
            w0_ = casadi::DM::zeros(Ndx + Nu, 1);
            general_lbw_ = -casadi::DM::inf(Ndx + Nu, 1);
            general_ubw_ = casadi::DM::inf(Ndx + Nu, 1);

            /*Transform initial guess for x to an initial guess for dx, using f_diff, the inverse of f_int*/
            casadi::MX xkg_sym = casadi::MX::sym("xkg", st_m_->nx, 1);
            casadi::MX xckg_sym = casadi::MX::sym("Xckg", st_m_->nx * dX_poly_.d, 1);
            if (!Wdata->initial_guess.is_null())
            {
                casadi::DM xg = Wdata->initial_guess.map(knot_num_ + 1, "serial")(knot_times_).at(0);
                casadi::Function dxg_func = casadi::Function("xg_fun", casadi::MXVector{xkg_sym}, casadi::MXVector{Fdiff_(casadi::MXVector{x0_global_, xkg_sym, 1.0}).at(0)})
                                                .map(knot_num_ + 1, "serial");
                w0_(casadi::Slice(0, Ndxknot)) = casadi::DM::reshape(dxg_func(casadi::DMVector{xg}).at(0), Ndxknot, 1);
                /*The transformation of xc to dxc is a slightly less trivial. While x_k = fint(x0_init, dx_k), for xc_k, we have xc_k = fint(x_k, dxc_k) which is equivalent to xc_k = fint(fint(x0_init, dx_k), dxc_k).
                Thus, dxc_k = fdiff(fint(x0_init, dx_k), xc_k)). This could be done with maps like above, but it is not necessary.*/
                casadi::DM xc_g = Wdata->initial_guess.map((dX_poly_.d) * knot_num_, "serial")(collocation_times_).at(0);
                for (casadi_int i = 0; i < knot_num_; ++i)
                {
                    casadi::DM xk = xg(casadi::Slice(i * st_m_->nx, (i + 1) * st_m_->nx));
                    casadi::DM xck = xc_g(casadi::Slice(i * st_m_->nx * dX_poly_.d, (i + 1) * st_m_->nx * dX_poly_.d));
                    for (casadi_int j = 0; j < dX_poly_.d; ++j)
                    {
                        w0_(casadi::Slice(Ndxknot + i * st_m_->ndx * dX_poly_.d + j * st_m_->ndx, Ndxknot + i * st_m_->ndx * dX_poly_.d + (j + 1) * st_m_->ndx)) =
                            reshape(Fdiff_(casadi::DMVector{xk, xck(casadi::Slice(j * st_m_->nx, (j + 1) * st_m_->nx)), h_}).at(0), st_m_->ndx, 1);
                    }
                }
            }

            if (!Wdata->lower_bound.is_null() && !Wdata->upper_bound.is_null())
            {
                // casadi::SX t = casadi::SX::sym("t");
                // casadi::Function lbdx = casadi::Function("lbdx", casadi::SXVector{t}, casadi::SXVector{Fdiff(casadi::SXVector{x0_global, Wdata->lower_bound(t).at(0), 1.0}).at(0)});
                // casadi::Function ubdx = casadi::Function("ubdx", casadi::SXVector{t}, casadi::SXVector{Fdiff(casadi::SXVector{x0_global, Wdata->upper_bound(t).at(0), 1.0}).at(0)});

                // general_lbw(casadi::Slice(0, Ndxknot)) = casadi::DM::reshape(lbdx.map(knot_num + 1, "serial")(knot_times).at(0), Ndxknot, 1);
                // general_ubw(casadi::Slice(0, Ndxknot)) = casadi::DM::reshape(ubdx.map(knot_num + 1, "serial")(knot_times).at(0), Ndxknot, 1);
                // general_lbw(casadi::Slice(Ndxknot, Ndx)) = casadi::DM::reshape(lbdx.map((dX_poly.d) * knot_num, "serial")(collocation_times).at(0), Ndxcol, 1);
                // general_ubw(casadi::Slice(Ndxknot, Ndx)) = casadi::DM::reshape(ubdx.map((dX_poly.d) * knot_num, "serial")(collocation_times).at(0), Ndxcol, 1);
                
                general_lbw_(casadi::Slice(0, Ndxknot)) = casadi::DM::reshape(Wdata->lower_bound.map(knot_num_ + 1, "serial")(knot_times_).at(0), Ndxknot, 1);
                general_ubw_(casadi::Slice(0, Ndxknot)) = casadi::DM::reshape(Wdata->upper_bound.map(knot_num_ + 1, "serial")(knot_times_).at(0), Ndxknot, 1);
                general_lbw_(casadi::Slice(Ndxknot, Ndx)) = casadi::DM::reshape(Wdata->lower_bound.map((dX_poly_.d) * knot_num_, "serial")(collocation_times_).at(0), Ndxcol, 1);
                general_ubw_(casadi::Slice(Ndxknot, Ndx)) = casadi::DM::reshape(Wdata->upper_bound.map((dX_poly_.d) * knot_num_, "serial")(collocation_times_).at(0), Ndxcol, 1);
            }
            if (!Wdata->initial_guess.is_null())
            {
                w0_(casadi::Slice(Ndx, Ndx + Nuknot)) = casadi::DM::reshape(Wdata->initial_guess.map(knot_num_ + 1, "serial")(u_knot_times_).at(1), Nuknot, 1);
                w0_(casadi::Slice(Ndx + Nuknot, Ndx + Nu)) = casadi::DM::reshape(Wdata->initial_guess.map((U_poly_.d) * knot_num_, "serial")(u_collocation_times_).at(1), Nucol, 1);
            }
            if (!Wdata->lower_bound.is_null() && !Wdata->upper_bound.is_null())
            {
                general_lbw_(casadi::Slice(Ndx, Ndx + Nuknot)) = casadi::DM::reshape(Wdata->lower_bound.map(knot_num_ + 1, "serial")(u_knot_times_).at(1), Nuknot, 1);
                general_ubw_(casadi::Slice(Ndx, Ndx + Nuknot)) = casadi::DM::reshape(Wdata->upper_bound.map(knot_num_ + 1, "serial")(u_knot_times_).at(1), Nuknot, 1);
                general_lbw_(casadi::Slice(Ndx + Nuknot, Ndx + Nu)) = casadi::DM::reshape(Wdata->lower_bound.map((U_poly_.d) * knot_num_, "serial")(u_collocation_times_).at(1), Nucol, 1);
                general_ubw_(casadi::Slice(Ndx + Nuknot, Ndx + Nu)) = casadi::DM::reshape(Wdata->upper_bound.map((U_poly_.d) * knot_num_, "serial")(u_collocation_times_).at(1), Nucol, 1);
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

            casadi::MX xs = processVector(X0_var_vec_);
            casadi::MX us = processVector(U0_var_vec_);
            casadi::MX dxs = processVector(dX0_var_vec_);
            casadi::MX dxcs = horzcat(dXc_var_vec_);
            casadi::MX ucs = horzcat(Uc_var_vec_);
            casadi::MX xs_offset = processOffsetVector(X0_var_vec_);
            casadi::MX dxs_offset = processOffsetVector(dX0_var_vec_);
            casadi::MX us_offset = processOffsetVector(U0_var_vec_);

            casadi::MXVector solmap_result = sol_map_func_(casadi::MXVector{xs, dxcs, dxs, us, ucs});
            casadi::MX all_xs = solmap_result.at(0);
            casadi::MX all_us = solmap_result.at(1);

            /*This section cannot get much faster, it is bounded by the time to evaluate the constraint*/
            casadi::MX col_con_mat = collocation_constraint_map_(casadi::MXVector{xs, dxcs, dxs, us, ucs}).at(0);
            casadi::MX xf_con_mat = xf_constraint_map_(casadi::MXVector{xs, dxcs, dxs, us, ucs}).at(0);
            casadi::MX uf_con_mat = uf_constraint_map_(casadi::MXVector{xs, dxcs, dxs, us, ucs}).at(0);
            dxs_offset = reshape(dxs_offset, dxs_offset.size1() * dxs_offset.size2(), 1);
            us_offset = reshape(us_offset, us_offset.size1() * us_offset.size2(), 1);

            result.push_back(reshape(col_con_mat, col_con_mat.size1() * col_con_mat.size2(), 1));
            result.push_back(reshape(xf_con_mat, xf_con_mat.size1() * xf_con_mat.size2(), 1) -
                             dxs_offset);
            result.push_back(reshape(uf_con_mat, uf_con_mat.size1() * uf_con_mat.size2(), 1) -
                             us_offset);

            for (size_t i = 0; i < general_constraint_maps_.size(); ++i)
            {
                casadi::MX g_con_mat = general_constraint_maps_[i](casadi::MXVector{xs, dxcs, dxs, us, ucs}).at(0);
                result.push_back(reshape(g_con_mat, g_con_mat.size1() * g_con_mat.size2(), 1));
            }

            casadi::MX cost = q_cost_fold_(casadi::MXVector{J0, xs, dxcs, dxs, us, ucs}).at(0);
            J0 = cost;
            /*where g of this segment starts*/
            size_t g_size = g.size();
            /*Use move to avoid copying the vectors. Reserve space for g in advance outside of PseudospectralSegment.*/
            g.insert(g.end(), make_move_iterator(result.begin()), make_move_iterator(result.end()));
            g_range_ = tuple_size_t(g_size, g.size());

            /*where w of this segment starts*/
            size_t w_it_size = w.size();
            size_t w_size = 0;
            /*Accumulate the size of each symbolic element in w*/
            if (w_it_size != 0)
                w_size = size_t(accumulate(w.begin(), w.end(), 0.0, [](int sum, const casadi::MX &item)
                                           { return sum + item.size1() * item.size2(); }));

            /*Use move to avoid copying the vectors. Reserve space for w in advance outside of PseudospectralSegment.*/
            w.insert(w.end(), make_move_iterator(dX0_var_vec_.begin()), make_move_iterator(dX0_var_vec_.end()));
            w.insert(w.end(), make_move_iterator(dXc_var_vec_.begin()), make_move_iterator(dXc_var_vec_.end()));
            w.insert(w.end(), make_move_iterator(U0_var_vec_.begin()), make_move_iterator(U0_var_vec_.end()));
            w.insert(w.end(), make_move_iterator(Uc_var_vec_.begin()), make_move_iterator(Uc_var_vec_.end()));

            w_range_ = tuple_size_t(w_size, accumulate(w.begin(), w.end(), 0.0, [](int sum, const casadi::MX &item)
                                                      { return sum + item.size1() * item.size2(); }));
            get_sol_func_ = casadi::Function("func",
                                            casadi::MXVector({vertcat(w)}),
                                            casadi::MXVector({all_xs, all_us}));
        }

        casadi::MXVector PseudospectralSegment::extractSolution(casadi::MX &w) const
        {
            return get_sol_func_(casadi::MXVector{w});
        }

        casadi::MX PseudospectralSegment::getInitialStateDeviant() const
        {
            return dX0_var_vec_.front();
        }

        casadi::MX PseudospectralSegment::getInitialState() const
        {
            return X0_var_vec_.front();
        }

        casadi::MX PseudospectralSegment::getFinalStateDeviant() const
        {
            return dX0_var_vec_.back();
        }

        casadi::MX PseudospectralSegment::getFinalState() const
        {
            return X0_var_vec_.back();
        }

        casadi::DM PseudospectralSegment::getSegmentTimes() const
        {
            return segment_times_;
        }

        casadi::DM PseudospectralSegment::getKnotTimes() const
        {
            return knot_times_;
        }

        casadi::DM PseudospectralSegment::getCollocationTimes() const
        {
            return collocation_times_;
        }

        casadi::DM PseudospectralSegment::getUSegmentTimes() const
        {
            return u_segment_times_;
        }

        casadi::DM PseudospectralSegment::getUKnotTimes() const
        {
            return u_knot_times_;
        }

        casadi::DM PseudospectralSegment::getUCollocationTimes() const
        {
            return u_collocation_times_;
        }

        void PseudospectralSegment::fill_lbw_ubw(std::vector<double> &lbw, std::vector<double> &ubw)
        {
            /*where lb/ub of this segment starts*/
            auto bw_size = lbw.size();
            std::vector<double> element_access1 = general_lbw_.get_elements();
            std::vector<double> element_access2 = general_ubw_.get_elements();

            lbw.insert(lbw.end(), element_access1.begin(), element_access1.end());
            ubw.insert(ubw.end(), element_access2.begin(), element_access2.end());
            lbw_ubw_range_ = tuple_size_t(bw_size, lbw.size());
        }

        void PseudospectralSegment::fill_lbg_ubg(std::vector<double> &lbg, std::vector<double> &ubg)
        {
            /*where lb/ub of this segment starts*/
            auto bg_size = lbg.size();
            std::vector<double> element_access1 = general_lbg_.get_elements();
            std::vector<double> element_access2 = general_ubg_.get_elements();

            lbg.insert(lbg.end(), element_access1.begin(), element_access1.end());
            ubg.insert(ubg.end(), element_access2.begin(), element_access2.end());

            lbg_ubg_range_ = tuple_size_t(bg_size, lbg.size());
        }

        void PseudospectralSegment::fill_w0(std::vector<double> &all_w0) const
        {
            std::vector<double> element_access1 = w0_.get_elements();
            all_w0.insert(all_w0.end(), element_access1.begin(), element_access1.end());
        }

        tuple_size_t PseudospectralSegment::get_range_idx_decision_variables() const
        {
            return w_range_;
        }

        tuple_size_t PseudospectralSegment::get_range_idx_constraint_expressions() const
        {
            return g_range_;
        }

        tuple_size_t PseudospectralSegment::get_range_idx_constraint_bounds() const
        {
            return lbg_ubg_range_;
        }

        tuple_size_t PseudospectralSegment::get_range_idx_decision_bounds() const
        {
            return lbw_ubw_range_;
        }

        tuple_size_t PseudospectralSegment::get_range_idx_time() const
        {
            return t_range_;
        }
    }
}