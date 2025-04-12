#ifndef __galileo_predictive_solvers_solver_base_hxx__
#define __galileo_predictive_solvers_solver_base_hxx__

#include "galileo/predictive/solvers/solver-base.hpp"

namespace galileo
{

    template <typename Derived>
    void SolverBase<Derived>::resizeData()
    {
        std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();
        for (std::size_t t = 0; t < num_segments; ++t)
        {
            std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
            const std::size_t nu = ocp_.trajectory_.phase_models_[phase_index].nu();
            const std::size_t ng = ocp_.trajectory_.phase_models_[phase_index].ng();
            us_[t].conservativeResize(nu);
            g_adj_[t].conservativeResize(ng);
        }
    }

    template <typename Derived>
    SolverBase<Derived>::NumScalar SolverBase<Derived>::computeDynamicFeasibility()
    {
        tmp_feas_ = 0.;
        const std::size_t num_phases = ocp_->get_num_phases();
        const Eigen::VectorXd &x0 = ocp_->get_x0();

        ocp_.trajectory_.phase_models_[0].stateDiff(xs_[0], x0, fs_[0]);
        std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();

        for (std::size_t t = 0; t < num_segments; ++t)
        {
            std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
            std::size_t segment_index = t - ocp_.trajectory_.phase_offsets_[phase_index];
            auto &phase_model = ocp_.trajectory_.phase_models_[phase_index];
            auto &phase_data = ocp_.trajectory_.phase_datas_[phase_index];
            phase_model.stateDiff(xs_[t + 1], phase_data.get_xnext_from_segment(segment_index), fs_[t + 1]);
        }

        switch (feasnorm_)
        {
        case LInf:
            tmp_feas_ = std::max(tmp_feas_, fs_[0].lpNorm<Eigen::Infinity>());
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                tmp_feas_ = std::max(tmp_feas_, fs_[t + 1].lpNorm<Eigen::Infinity>());
            }
            break;
        case L1:
            tmp_feas_ = fs_[0].lpNorm<1>();
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                tmp_feas_ += fs_[t + 1].lpNorm<1>();
            }
            break;
        }
        return tmp_feas_;
    }

    template <typename Derived>
    SolverBase<Derived>::NumScalar SolverBase<Derived>::computeEqualityFeasibility()
    {
        tmp_feas_ = 0.;
        const std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();
        switch (feasnorm_)
        {
        case LInf:
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                auto &phase_model = ocp_.trajectory_.phase_models_[phase_index];
                auto &phase_data = ocp_.trajectory_.phase_datas_[phase_index];
                if (phase_model.nh() > 0)
                {
                    tmp_feas_ =
                        std::max(tmp_feas_, phase_data.get_H_from_segment(t).lpNorm<Eigen::Infinity>());
                }
            }
            break;
        case L1:
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                auto &phase_model = ocp_.trajectory_.phase_models_[phase_index];
                auto &phase_data = ocp_.trajectory_.phase_datas_[phase_index];
                if (phase_model.nh() > 0)
                {
                    tmp_feas_ += phase_data.get_H_from_segment(t).lpNorm<1>();
                }
            }
            break;
        }
        return tmp_feas_;
    }

    template <typename Derived>
    SolverBase<Derived>::NumScalar SolverBase<Derived>::computeInequalityFeasibility()
    {
        tmp_feas_ = 0.;
        const std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();
        switch (feasnorm_)
        {
        case LInf:
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                auto &phase_model = ocp_.trajectory_.phase_models_[phase_index];
                auto &phase_data = ocp_.trajectory_.phase_datas_[phase_index];
                if (phase_model.ng() > 0)
                {
                    g_adj_[t] = phase_data.get_G_from_segment(t)
                                    .cwiseMax(phase_model.get_g_lb_from_segment(t))
                                    .cwiseMin(phase_model.get_g_ub_from_segment(t));
                    tmp_feas_ = std::max(
                        tmp_feas_, (phase_data.get_G_from_segment(t) - g_adj_[t]).lpNorm<Eigen::Infinity>());
                }
            }
            break;
        }
        return tmp_feas_;
    }

    template <typename Derived>
    void SolverBase<Derived>::setCandidate(const std::vector<VectorXn> &xs_warm,
                                           const std::vector<VectorXn> &us_warm,
                                           bool is_feasible = false)
    {
        const std::size_t num_segments = ocp_.trajectory_.phase_offsets_.back();

        if (xs_warm.size() == 0)
        {
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
                xs_[t] = ocp_.trajectory_.phase_models_[phase_index].stateZero();
            }
            xs_.back() = ocp_.trajectory_.phase_models_.back().stateZero();
        }
        else
        {
            std::copy(xs_warm.begin(), xs_warm.end(), xs_.begin());
        }

        if (us_warm.size() == 0)
        {
            for (std::size_t t = 0; t < num_segments; ++t)
            {
                std::size_t phase_index = ocp_.trajectory_.find_phase_index(t);
                const std::size_t nu = ocp_.trajectory_.phase_models_[phase_index].nu();
                us_[t] = Eigen::VectorXd::Zero(nu);
            }
        }
        else
        {
            std::copy(us_warm.begin(), us_warm.end(), us_.begin());
        }
        is_feasible_ = is_feasible;
    }

} // namespace galileo

#endif // __galileo_predictive_solvers_solver_base_hxx__