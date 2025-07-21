#ifndef __galileo_predictive_segments_segment_erk_model_base_hpp__
#define __galileo_predictive_segments_segment_erk_model_base_hpp__

#include "galileo/predictive/phases/phase-spec.hpp"
#include "galileo/predictive/segments/segment-erk-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class SegmentERKModelBase
        : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using NumScalar = typename PS::NumScalar;
        using State_t = typename PS::State_t;

        template <typename StateVectorType, typename ControlParamVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            this->derived().calc(data, x.derived(), w.derived());
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x.derived());
        }

        template <typename StateVectorType, typename ControlParamVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            this->derived().calcDiff(data, x.derived(), w.derived());
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x.derived());
        }

        template <typename StateVectorType, typename ControlParamVectorType>
        void quasiStatic(Data_t &data,
                         const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlParamVectorType> &w,
                         const int maxiter,
                         const NumScalar &tol) const
        {
            this->derived().quasiStatic(data, x.derived(), w.derived(), maxiter, tol);
        }

        Data_t createData()
        {
            return this->derived().createData();
        }

        const PS &get_ps() const
        {
            return ps_.get();
        }

        const State_t &get_state() const
        {
            return state_.get();
        }

    protected:
        inline SegmentERKModelBase(const PS &ps)
            : ps_(ps), state_(ps.get_state())
        {
        }

        inline SegmentERKModelBase(const SegmentERKModelBase &clone)
            : ps_(clone.ps_), state_(clone.state_)
        {
        }

        inline SegmentERKModelBase &operator=(const SegmentERKModelBase &clone)
        {
            ps_ = clone.ps_;
            state_ = clone.state_;
            return *this;
        }

        std::reference_wrapper<const PS> ps_;
        std::reference_wrapper<const State_t> state_;

    }; // class SegmentERKModelBase

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_model_base_hpp__
