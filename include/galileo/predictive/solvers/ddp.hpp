#ifndef __galileo_predictive_solvers_ddp_hpp__
#define __galileo_predictive_solvers_ddp_hpp__

#include "galileo/predictive/solvers/fwd.hpp"
#include "galileo/predictive/solvers/solver-base.hpp"

#include <cmath>
#include <iostream>
#include <string>

namespace galileo
{

    template <typename _VarScalar, typename _NumScalar, int _Options, FeasibilityNormOptions _FeasibilityNorm,
              template <typename V, typename N, int O> class PhaseCollectionTpl>
    class SolverDDP; // forward declaration

    template <typename _VarScalar, typename _NumScalar, int _Options, FeasibilityNormOptions _FeasibilityNorm,
              template <typename V, typename N, int O> class PhaseCollectionTpl>
    struct traits<SolverDDP<_VarScalar, _NumScalar, _Options, _FeasibilityNorm, PhaseCollectionTpl>>
    {
        using SolverDerived = SolverDDP<_VarScalar, _NumScalar, _Options, _FeasibilityNorm, PhaseCollectionTpl>;
        using VarScalar = _VarScalar;
        using NumScalar = _NumScalar;
        static constexpr int Options = _Options;
        static constexpr FeasibilityNormOptions FeasibilityNorm = _FeasibilityNorm;
        using OptimalControlProblem_t = OptimalControlProblem<VarScalar, NumScalar, Options, PhaseCollectionTpl>;
        using VectorXv = Eigen::GMatrix<VarScalar, Eigen::Dynamic, 1>;
        using VectorXn = Eigen::GMatrix<NumScalar, Eigen::Dynamic, 1>;
    }; // struct traits

    template <typename _VarScalar, typename _NumScalar, int _Options, FeasibilityNormOptions _FeasibilityNorm,
              template <typename, typename, int> class PhaseCollectionTpl>
    class SolverDDP : SolverBase<SolverDDP<_VarScalar, _NumScalar, _Options, _FeasibilityNorm, PhaseCollectionTpl>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using SolverDerived = SolverDDP<_VarScalar, _NumScalar, _Options, _FeasibilityNorm, PhaseCollectionTpl>;
        GALILEO_SOLVER_BASIC_TYPEDEF(SolverDerived);
        GALILEO_SOLVER_TYPEDEF(SolverDerived);

        using MatrixXv = Eigen::GMatrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic>;
        using MatrixXn = Eigen::GMatrix<NumScalar, Eigen::Dynamic, Eigen::Dynamic>;
        using MatrixXvRowMajor = Eigen::GMatrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using MatrixXnRowMajor = Eigen::GMatrix<NumScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using Vector2n = Eigen::GMatrix<NumScalar, 2, 1>;

        SolverDDP(OptimalControlProblem_t ocp)
            : ocp_(ocp),
              reg_incfactor_(10.),
              reg_decfactor_(10.),
              reg_min_(1e-9),
              reg_max_(1e9),
              cost_try_(0.),
              th_grad_(1e-12),
              th_stepdec_(0.5),
              th_stepinc_(0.01)
        {
            allocateData();

            const std::size_t n_alphas = 10;
            alphas_.resize(n_alphas);
            for (std::size_t n = 0; n < n_alphas; ++n)
            {
                alphas_[n] = 1. / pow(2., static_cast<double>(n));
            }
            if (th_stepinc_ < alphas_[n_alphas - 1])
            {
                th_stepinc_ = alphas_[n_alphas - 1];
                std::cerr << "Warning: th_stepinc has higher value than lowest alpha "
                             "value, set to "
                          << std::to_string(alphas_[n_alphas - 1]) << std::endl;
            }
        }

        bool solve(
            const std::vector<VectorXn> &init_xs,
            const std::vector<VectorXn> &init_us,
            const std::size_t maxiter = 100,
            const bool is_feasible = false,
            const NumScalar init_reg = NAN)
        {
            if (ocp_.is_updated())
            {
                resizeData();
            }
            xs_try_[0] = ocp_.get_x0(); // it is needed in case that init_xs[0] is infeasible
            setCandidate(init_xs, init_us, is_feasible);

            if (std::isnan(init_reg))
            {
                preg_ = reg_min_;
                dreg_ = reg_min_;
            }
            else
            {
                preg_ = init_reg;
                dreg_ = init_reg;
            }
            was_feasible_ = false;

            bool recalcDiff = true;
            for (iter_ = 0; iter_ < maxiter; ++iter_)
            {
                while (true)
                {
                    try
                    {
                        computeDirection(recalcDiff);
                    }
                    catch (std::exception &e)
                    {
                        recalcDiff = false;
                        increaseRegularization();
                        if (preg_ == reg_max_)
                        {
                            return false;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    break;
                }
                expectedImprovement();

                // We need to recalculate the derivatives when the step length passes
                recalcDiff = false;
                for (std::vector<NumScalar>::const_iterator it = alphas_.begin();
                     it != alphas_.end(); ++it)
                {
                    steplength_ = *it;

                    try
                    {
                        dV_ = tryStep(steplength_);
                    }
                    catch (std::exception &e)
                    {
                        continue;
                    }
                    dVexp_ = steplength_ * (d_[0] + 0.5 * steplength_ * d_[1]);

                    if (dVexp_ >= 0)
                    { // descend direction
                        if (std::abs(d_[0]) < th_grad_ || !is_feasible_ ||
                            dV_ > th_acceptstep_ * dVexp_)
                        {
                            was_feasible_ = is_feasible_;
                            setCandidate(xs_try_, us_try_, true);
                            cost_ = cost_try_;
                            recalcDiff = true;
                            break;
                        }
                    }
                }

                if (steplength_ > th_stepdec_)
                {
                    decreaseRegularization();
                }
                if (steplength_ <= th_stepinc_)
                {
                    increaseRegularization();
                    if (preg_ == reg_max_)
                    {
                        return false;
                    }
                }
                stoppingCriteria();

                if (was_feasible_ && stop_ < th_stop_)
                {
                    return true;
                }
            }
            return false;
        }

        void computeDirection(const bool recalcDiff = true)
        {
            if (recalcDiff)
            {
                calcDiff();
            }
            backwardPass();
        }

        NumScalar tryStep(const NumScalar steplength = 1)
        {
            forwardPass(steplength);
            return cost_ - cost_try_;
        }

        NumScalar stoppingCriteria()
        {
            stop_ = std::abs(d_[0] + 0.5 * d_[1]);
            return stop_;
        }

        const Vector2n &expectedImprovement()
        {
            d_.fill(0);
            const std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
                const std::size_t nu = ocp_.trajectory_.phase_models_[phase_index].nu();
                if (nu != 0)
                {
                    d_[0] += Qu_[t].dot(k_[t]);
                    d_[1] -= k_[t].dot(Quuk_[t]);
                }
            }
            return d_;
        }

        void resizeData()
        {
            SolverBase::resizeData();

            const std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
                const std::size_t ndx = ocp_.trajectory_.phase_models_[phase_index].ndx();
                const std::size_t nu = ocp_.trajectory_.phase_models_[phase_index].nu();
                Qxu_[t].conservativeResize(ndx, nu);
                Quu_[t].conservativeResize(nu, nu);
                Qu_[t].conservativeResize(nu);
                K_[t].conservativeResize(nu, ndx);
                k_[t].conservativeResize(nu);
                us_try_[t].conservativeResize(nu);
                FuTVxx_p_[t].conservativeResize(nu, ndx);
                Quuk_[t].conservativeResize(nu);
                if (nu != 0)
                {
                    FuTVxx_p_[t].setZero();
                }
            }
        }

        NumScalar calcDiff()
        {
            if (iter_ == 0)
            {
                ocp_.calc(xs_, us_);
            }
            cost_ = ocp_.calcDiff(xs_, us_);

            ffeas_ = computeDynamicFeasibility();
            gfeas_ = computeInequalityFeasibility();
            hfeas_ = computeEqualityFeasibility();
            return cost_;
        }

        void backwardPass()
        {
            if (!std::isnan(preg_))
            {
                Vxx_.back().diagonal().array() += preg_;
            }

            if (!is_feasible_)
            {
                Vx_.back().noalias() += Vxx_.back() * fs_.back();
            }

            const std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();
            for (std::size_t t = num_segments; t > 0; --t)
            {
                // Compute the linear-quadratic approximation of the control Hamiltonian
                // function
                computeActionValueFunction(t);

                // Compute the feedforward and feedback gains
                computeGains(t);

                // Compute the linear-quadratic approximation of the Value function
                computeValueFunction(t);
            }
        }

        void forwardPass(const NumScalar stepLength)
        {
            cost_try_ = 0.;
            const std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
                auto &phase_model = ocp_.trajectory_.phase_models_[phase_index];
                auto &phase_data = ocp_.trajectory_.phase_datas_[phase_index];

                phase_model.stateDiff(xs_[t], xs_try_[t], dx_[t]);
                if (phase_model.nu() != 0)
                {
                    us_try_[t].noalias() = us_[t];
                    us_try_[t].noalias() -= k_[t] * steplength;
                    us_try_[t].noalias() -= K_[t] * dx_[t];
                    phase_model.segmentCalc(t, xs_try_[t], us_try_[t]);
                }
                else
                {
                    phase_model.segmentCalc(t, xs_try_[t]);
                }
                xs_try_[t + 1] = phase_data.get_xnext_from_segment(t);
                cost_try_ += phase_data.get_cost_from_segment(t);
            }
        }

        void computeActionValueFunction(
            const std::size_t t)
        {
            std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
            const std::size_t nu = ocp_.trajectory_.phase_models_[phase_index].nu();
            const Eigen::MatrixXd &Vxx_p = Vxx_[t + 1];
            const Eigen::VectorXd &Vx_p = Vx_[t + 1];

            auto &phase_data = ocp_.trajectory_.phase_datas_[phase_index];

            FxTVxx_p_.noalias() = phase_data.get_Fx_from_segment(t).transpose() * Vxx_p;
            Qx_[t] = phase_data.get_Lx_from_segment(t);
            Qx_[t].noalias() += phase_data.get_Fx_from_segment(t).transpose() * Vx_p;
            Qxx_[t] = phase_data.get_Lxx_from_segment(t);
            Qxx_[t].noalias() += FxTVxx_p_ * phase_data.get_Fx_from_segment(t);
            if (nu != 0)
            {
                FuTVxx_p_[t].noalias() = phase_data.get_Fu_from_segment(t).transpose() * Vxx_p;
                Qu_[t] = phase_data.get_Lu_from_segment(t);
                Qu_[t].noalias() += phase_data.get_Fu_from_segment(t).transpose() * Vx_p;
                Quu_[t] = phase_data.get_Luu_from_segment(t);
                Quu_[t].noalias() += FuTVxx_p_[t] * phase_data.get_Fu_from_segment(t);
                Qxu_[t] = phase_data.get_Lxu_from_segment(t);
                Qxu_[t].noalias() += FxTVxx_p_ * phase_data.get_Fu_from_segment(t);
                if (!std::isnan(preg_))
                {
                    Quu_[t].diagonal().array() += preg_;
                }
            }
        }

        void computeValueFunction(const std::size_t t)
        {
            std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
            const std::size_t nu = ocp_.trajectory_.phase_models_[phase_index].nu();
            Vx_[t] = Qx_[t];
            Vxx_[t] = Qxx_[t];
            if (nu != 0)
            {
                Quuk_[t].noalias() = Quu_[t] * k_[t];
                Vx_[t].noalias() -= K_[t].transpose() * Qu_[t];
                Vxx_[t].noalias() -= Qxu_[t] * K_[t];
            }
            Vxx_tmp_ = 0.5 * (Vxx_[t] + Vxx_[t].transpose());
            Vxx_[t] = Vxx_tmp_;

            if (!std::isnan(preg_))
            {
                Vxx_[t].diagonal().array() += preg_;
            }

            // Compute and store the Vx gradient at end of the interval (rollout state)
            if (!is_feasible_)
            {
                Vx_[t].noalias() += Vxx_[t] * fs_[t];
            }
        }

        void computeGains(const std::size_t t)
        {
            std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
            const std::size_t nu = ocp_.trajectory_.phase_models_[phase_index].nu();
            if (nu > 0)
            {
                Quu_llt_[t].compute(Quu_[t]);
                const Eigen::ComputationInfo &info = Quu_llt_[t].info();
                if (info != Eigen::Success)
                {
                    throw std::runtime_error("backward_error");
                }
                K_[t] = Qxu_[t].transpose();

                Quu_llt_[t].solveInPlace(K_[t]);
                k_[t] = Qu_[t];
                Quu_llt_[t].solveInPlace(k_[t]);
            }
        }

        void increaseRegularization()
        {
            preg_ *= reg_incfactor_;
            if (preg_ > reg_max_)
            {
                preg_ = reg_max_;
            }
            dreg_ = preg_;
        }

        void decreaseRegularization()
        {
            preg_ /= reg_decfactor_;
            if (preg_ < reg_min_)
            {
                preg_ = reg_min_;
            }
            dreg_ = preg_;
        }

        void allocateData()
        {
            const std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();
            Vxx_.resize(num_segments + 1);
            Vx_.resize(num_segments + 1);
            Qxx_.resize(num_segments);
            Qxu_.resize(num_segments);
            Quu_.resize(num_segments);
            Qx_.resize(num_segments);
            Qu_.resize(num_segments);
            K_.resize(num_segments);
            k_.resize(num_segments);

            xs_try_.resize(num_segments + 1);
            us_try_.resize(num_segments);
            dx_.resize(num_segments);

            FuTVxx_p_.resize(num_segments);
            Quu_llt_.resize(num_segments);
            Quuk_.resize(num_segments);

            for (std::size_t t = 0; t < num_segments; ++t)
            {
                std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
                const std::size_t ndx = ocp_.trajectory_.phase_models_[phase_index].ndx();
                const std::size_t nu = ocp_.trajectory_.phase_models_[phase_index].nu();
                Vxx_[t] = Eigen::MatrixXd::Zero(ndx, ndx);
                Vx_[t] = Eigen::VectorXd::Zero(ndx);
                Qxx_[t] = Eigen::MatrixXd::Zero(ndx, ndx);
                Qxu_[t] = Eigen::MatrixXd::Zero(ndx, nu);
                Quu_[t] = Eigen::MatrixXd::Zero(nu, nu);
                Qx_[t] = Eigen::VectorXd::Zero(ndx);
                Qu_[t] = Eigen::VectorXd::Zero(nu);
                K_[t] = MatrixXdRowMajor::Zero(nu, ndx);
                k_[t] = Eigen::VectorXd::Zero(nu);

                if (t == 0)
                {
                    xs_try_[t] = ocp_.trajectory_.phase_models_[phase_index].stateZero();
                }
                else
                {
                    xs_try_[t] = ocp_.trajectory_.phase_models_[phase_index].stateZero();
                }
                us_try_[t] = Eigen::VectorXd::Zero(nu);
                dx_[t] = Eigen::VectorXd::Zero(ndx);

                FuTVxx_p_[t] = MatrixXdRowMajor::Zero(nu, ndx);
                Quu_llt_[t] = Eigen::LLT<Eigen::MatrixXd>(nu);
                Quuk_[t] = Eigen::VectorXd(nu);
            }
            Vxx_.back() = Eigen::MatrixXd::Zero(ndx, ndx);
            Vxx_tmp_ = Eigen::MatrixXd::Zero(ndx, ndx);
            Vx_.back() = Eigen::VectorXd::Zero(ndx);
            xs_try_.back() = ocp_.trajectory_.phase_models_.back().stateZero();

            FxTVxx_p_ = MatrixXdRowMajor::Zero(ndx, ndx);
            fTVxx_p_ = Eigen::VectorXd::Zero(ndx);
        }

    protected:
        OptimalControlProblem_t ocp_;

        NumScalar merit_;
        NumScalar stop_;

        Vector2n d_;
        NumScalar dV_;
        NumScalar dPhi_;
        NumScalar dVexp_;
        NumScalar dPhiexp_;
        NumScalar dfeas_;

        NumScalar reg_incfactor_;
        NumScalar reg_decfactor_;

        NumScalar reg_min_;
        NumScalar reg_max_;

        NumScalar cost_try_;
        std::vector<VectorXn> xs_try_;
        std::vector<VectorXn> us_try_;
        std::vector<VectorXn> dx_;

        std::vector<MatrixXn> Vxx_;
        MatrixXn Vxx_tmp_;
        std::vector<VectorXn> Vx_;
        std::vector<MatrixXn> Qxx_;
        std::vector<MatrixXn> Qxu_;
        std::vector<MatrixXn> Quu_;
        std::vector<VectorXn> Qx_;
        std::vector<VectorXn> Qu_;
        std::vector<MatrixXnRowMajor> K_;
        std::vector<VectorXn> k_;

        VectorXn xnext_;
        MatrixXnRowMajor FxTVxx_p_;

        std::vector<MatrixXnRowMajor> FuTVxx_p_;

        VectorXn fTVxx_p_;

        std::vector<Eigen::LLT<MatrixXn>> Quu_llt_;
        std::vector<VectorXn> Quuk_;

        std::vector<NumScalar> alphas_;
        NumScalar th_grad_;

        NumScalar th_stepdec_;
        NumScalar th_stepinc_;

    }; // class SolverDDP

} // namespace galileo

// #include "galileo/predictive/solvers/ddp.hxx"

#endif // __galileo_predictive_solvers_ddp_hpp__
