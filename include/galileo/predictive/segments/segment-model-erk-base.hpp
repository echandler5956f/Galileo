#ifndef __galileo_predictive_segments_segment_model_erk_base_hpp__
#define __galileo_predictive_segments_segment_model_erk_base_hpp__

#include "galileo/predictive/segments/segment-model-base.hpp"

namespace galileo
{
    namespace predictive
    {

        template <typename Derived>
        class SegmentModelERKBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using NodeDerived = typename traits<Derived>::NodeDerived;
            using SegmentDerived = typename traits<Derived>::SegmentDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_BASIC_TYPEDEF(SegmentDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_CONSTANTS(SegmentDerived);

            GALILEO_NODE_MODEL_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_MODEL_ERK_TYPEDEF(SegmentDerived);

            template <typename StateVectorType, typename ControlMatrixType>
            void calc(SegmentDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calc(data, x.derived(), us.derived());
            }

            template <typename StateVectorType, typename ControlMatrixType>
            void calcDiff(SegmentDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calcDiff(data, x.derived(), us.derived());
            }

            template <typename StateVectorType, typename ControlMatrixType>
            void quasiStatic(SegmentDataDerived &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlMatrixType> &us,
                             const std::size_t maxiter, const Scalar tol) const
            {
                derived().quasiStatic(data, x.derived(), us.derived(), maxiter, tol);
            }

        protected:
            inline SegmentModelERKBase()
            {
            }

            inline SegmentModelERKBase(const SegmentModelERKBase &clone)
            {
                *this = clone;
            }

            inline SegmentModelERKBase &operator=(const SegmentModelERKBase &clone)
            {
                return *this;
            }

        }; // class SegmentModelERKBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_model_erk_base_hpp__
