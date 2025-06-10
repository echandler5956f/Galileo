#ifndef __galileo_predictive_segments_segment_erk_model_base_hpp__
#define __galileo_predictive_segments_segment_erk_model_base_hpp__

#include "galileo/predictive/segments/segment-erk-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class SegmentERKModelBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

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
                         const std::size_t maxiter,
                         const typename PS::NumScalar &tol) const
        {
            this->derived().quasiStatic(data, x.derived(), w.derived(), maxiter, tol);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector)
        {
            return this->derived().createData(collector);
        }

    protected:
        inline SegmentERKModelBase()
        {
        }

        inline SegmentERKModelBase(const SegmentERKModelBase &clone)
        {
            *this = clone;
        }

        inline SegmentERKModelBase &operator=(const SegmentERKModelBase &clone)
        {
            return *this;
        }

    }; // class SegmentERKModelBase

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_model_base_hpp__
