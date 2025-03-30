#ifndef __galileo_predictive_segments_segment_model_erk_base_hpp__
#define __galileo_predictive_segments_segment_model_erk_base_hpp__

#include "galileo/predictive/segments/segment-base.hpp"
#include "galileo/predictive/nodes/node-base.hpp"

namespace galileo
{
    namespace predictive
    {

        enum class ERKType
        {
            Euler = 1,
            RK2 = 2,
            RK3 = 3,
            RK4 = 4
        };

        template <typename Derived, typename PhaseSpec>
        class SegmentModelERKBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            template <typename StateVectorType, typename ControlMatrixType>
            void calc(typename PS::SegmentData_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calc(data, x.derived(), us.derived());
            }

            template <typename StateVectorType, typename ControlMatrixType>
            void calcDiff(typename PS::SegmentData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calcDiff(data, x.derived(), us.derived());
            }

            template <typename StateVectorType, typename ControlMatrixType>
            void quasiStatic(typename PS::SegmentData_t &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlMatrixType> &us,
                             const std::size_t maxiter, const typename PS::NumScalar tol) const
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
