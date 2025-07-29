#ifndef __galileo_predictive_phases_phase_model_base_hpp__
#define __galileo_predictive_phases_phase_model_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class PhaseModelBase
        : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs,
                  const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            int i = 0;
            for (auto model_it = segments_.begin(), data_it = data.segments.begin();
                 model_it != segments_.end(); ++model_it, ++data_it, ++i)
            {
                model_it->calc(*data_it, col(xs, i), col(ws, i));
            }
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            int i = 0;
            for (auto model_it = segments_.begin(), data_it = data.segments.begin();
                 model_it != segments_.end(); ++model_it, ++data_it, ++i)
            {
                model_it->calcDiff(*data_it, col(xs, i), col(ws, i));
            }
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void quasiStatic(Data_t &data, const Eigen::MatrixBase<StateMatrixType> &xs,
                         Eigen::MatrixBase<ControlParamMatrixType> &ws,
                         const int maxiter, const NumScalar &tol) const
        {
            int i = 0;
            for (auto model_it = segments_.begin(), data_it = data.segments.begin();
                 model_it != segments_.end(); ++model_it, ++data_it, ++i)
            {
                model_it->quasiStatic(*data_it, col(xs, i), col(ws, i), maxiter, tol);
            }
        }

        template <bool IsForward, typename RightPhaseModelType, typename RightPhaseDataType, typename StateVectorType, typename ControlVectorType>
        void resetMapCalc(Data_t &data,
                          const RightPhaseModelType &right_model,
                          RightPhaseDataType &right_data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().resetMapCalc(data, right_model, right_data, x, u);
        }

        template <bool IsForward, typename RightPhaseModelType, typename RightPhaseDataType, typename StateVectorType>
        void resetMapCalc(Data_t &data,
                          const RightPhaseModelType &right_model,
                          RightPhaseDataType &right_data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().resetMapCalc(data, right_model, right_data, x);
        }

        template <bool IsForward, typename RightPhaseModelType, typename RightPhaseDataType, typename StateVectorType, typename ControlVectorType>
        void resetMapCalcDiff(Data_t &data,
                              const RightPhaseModelType &right_model,
                              RightPhaseDataType &right_data,
                              const Eigen::MatrixBase<StateVectorType> &x,
                              const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().resetMapCalcDiff(data, right_model, right_data, x, u);
        }

        template <bool IsForward, typename RightPhaseModelType, typename RightPhaseDataType, typename StateVectorType>
        void resetMapCalcDiff(Data_t &data,
                              const RightPhaseModelType &right_model,
                              RightPhaseDataType &right_data,
                              const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().resetMapCalcDiff(data, right_model, right_data, x);
        }

        Data_t createData() const
        {
            return Data_t(*this);
        }

        const PS &get_ps() const
        {
            return this->derived().get_ps_impl();
        }

        const PS &get_ps_impl() const
        {
            return ps_.get();
        }

        const std::vector<SegmentModel_t> &get_segments() const
        {
            return segments_;
        }

    protected:
        inline PhaseModelBase(const PS &ps)
            : ps_(ps)
        {
        }

        inline PhaseModelBase(const PhaseModelBase &clone)
            : ps_(clone.ps_)
        {
            *this = clone;
        }

        inline PhaseModelBase &operator=(const PhaseModelBase &clone)
        {
            ps_ = clone.ps_;
            return *this;
        }

        std::reference_wrapper<const PS> ps_;
        std::vector<SegmentModel_t> segments_;

    }; // class PhaseModelBase

} // namespace galileo

#endif // __galileo_predictive_phases_phase_model_base_hpp__
