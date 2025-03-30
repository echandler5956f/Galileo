#ifndef __galileo_predictive_segments_segment_model_erk_base_hpp__
#define __galileo_predictive_segments_segment_model_erk_base_hpp__

#include "galileo/predictive/segments/segment-base.hpp"
#include "galileo/predictive/nodes/node-base.hpp"

#define GALILEO_SEGMENT_ERK_BASIC_TYPEDEF(Segment)                                   \
    using SegmentModelERKDerived = typename traits<Segment>::SegmentModelERKDerived; \
    using SegmentDataERKDerived = typename traits<Segment>::SegmentDataERKDerived;

#define GALILEO_SEGMENT_ERK_CONSTANTS(Segment) \
    static constexpr int NStages = traits<Segment>::NStages;

#define GALILEO_SEGMENT_MODEL_ERK_TYPEDEF(Segment)                                       \
    using NodeModelDerived = typename traits<Segment>::NodeModelDerived;                 \
    using ControlParamModelDerived = typename traits<Segment>::ControlParamModelDerived; \
    using StageCoefficients_t = typename traits<Segment>::StageCoefficients_t;           \
    using Quadrature_t = typename traits<Segment>::Quadrature_t;                         \
    using Timings_t = typename traits<Segment>::Timings_t;

#define GALILEO_SEGMENT_DATA_ERK_TYPEDEF(Segment)                          \
    using NodeDataVector = typename traits<Segment>::NodeDataVector;       \
    using ControlDataVector = typename traits<Segment>::ControlDataVector; \
    using F_t = typename traits<Segment>::F_t;                             \
    using Fx_t = typename traits<Segment>::Fx_t;                           \
    using Fu_t = typename traits<Segment>::Fu_t;                           \
    using L_t = typename traits<Segment>::L_t;                             \
    using Lx_t = typename traits<Segment>::Lx_t;                           \
    using Lu_t = typename traits<Segment>::Lu_t;                           \
    using Lxx_t = typename traits<Segment>::Lxx_t;                         \
    using Lxu_t = typename traits<Segment>::Lxu_t;                         \
    using Luu_t = typename traits<Segment>::Luu_t;                         \
    using H_t = typename traits<Segment>::H_t;                             \
    using Hx_t = typename traits<Segment>::Hx_t;                           \
    using Hu_t = typename traits<Segment>::Hu_t;                           \
    using G_t = typename traits<Segment>::G_t;                             \
    using Gx_t = typename traits<Segment>::Gx_t;                           \
    using Gu_t = typename traits<Segment>::Gu_t;

namespace galileo
{
    namespace predictive
    {

        enum ERKType : int
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

            using NodeDerived = typename traits<Derived>::NodeDerived;
            using SegmentDerived = typename traits<Derived>::SegmentDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_ERK_BASIC_TYPEDEF(SegmentDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_ERK_CONSTANTS(SegmentDerived);

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
