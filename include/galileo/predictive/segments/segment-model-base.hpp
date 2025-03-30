#ifndef __galileo_predictive_segments_segment_model_base_hpp__
#define __galileo_predictive_segments_segment_model_base_hpp__

#include "galileo/predictive/segments/segment-base.hpp"
#include "galileo/predictive/nodes/node-base.hpp"

namespace galileo
{
    namespace predictive
    {

        template <typename Derived, typename PhaseSpec>
        class SegmentModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            template <typename StateMatrixType, typename ControlMatrixType>
            void calc(typename PS::SegmentData_t &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calc(data, xs.derived(), us.derived());
            }

            template <typename StateMatrixType, typename ControlMatrixType>
            void calcDiff(typename PS::SegmentData_t &data,
                          const Eigen::MatrixBase<StateMatrixType> &xs,
                          const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calcDiff(data, xs.derived(), us.derived());
            }

            template <typename StateMatrixType, typename ControlMatrixType>
            void quasiStatic(typename PS::SegmentData_t &data, const Eigen::MatrixBase<StateMatrixType> &xs,
                             Eigen::MatrixBase<ControlMatrixType> &us,
                             const std::size_t maxiter, const typename PS::NumScalar tol) const
            {
                derived().quasiStatic(data, xs.derived(), us.derived(), maxiter, tol);
            }

        protected:
            inline SegmentModelBase()
            {
            }

            inline SegmentModelBase(const SegmentModelBase &clone)
            {
                *this = clone;
            }

            inline SegmentModelBase &operator=(const SegmentModelBase &clone)
            {
                return *this;
            }

        }; // class SegmentModelBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_model_base_hpp__
