#ifndef __galileo_predictive_segments_segment_erk_model_base_hpp__
#define __galileo_predictive_segments_segment_erk_model_base_hpp__

#include "galileo/predictive/segments/segment-erk-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{
    namespace predictive
    {

        template <typename Derived, typename PhaseSpec>
        class SegmentERKModelBase : internal::CRTP<SegmentERKModelBase<Derived, PhaseSpec>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using SegmentERKDerived = typename traits<Derived>::SegmentERKDerived;
            using SegmentERKDataDerived = typename traits<SegmentERKDerived>::SegmentERKDataDerived;
            using SegmentERKModelDerived = typename traits<SegmentERKDerived>::SegmentERKModelDerived;

            template <typename StateVectorType, typename ControlParamVectorType>
            void calc(SegmentERKDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
            {
                derived().calc(data, x.derived(), w.derived());
            }

            template <typename StateVectorType>
            void calc(SegmentERKDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calc(data, x.derived());
            }

            template <typename StateVectorType, typename ControlParamVectorType>
            void calcDiff(SegmentERKDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlParamVectorType> &w) const
            {
                derived().calcDiff(data, x.derived(), w.derived());
            }

            template <typename StateVectorType>
            void calcDiff(SegmentERKDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calcDiff(data, x.derived());
            }

            template <typename StateVectorType, typename ControlParamVectorType>
            void quasiStatic(SegmentERKDataDerived &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlParamVectorType> &w,
                             const std::size_t maxiter,
                             const typename PS::NumScalar &tol) const
            {
                derived().quasiStatic(data, x.derived(), w.derived(), maxiter, tol);
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

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_model_base_hpp__
