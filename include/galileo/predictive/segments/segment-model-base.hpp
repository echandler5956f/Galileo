#ifndef __galileo_predictive_segments_segment_model_base_hpp__
#define __galileo_predictive_segments_segment_model_base_hpp__

#include "galileo/predictive/segments/segment-base.hpp"
#include "galileo/predictive/nodes/node-base.hpp"

#define GALILEO_SEGMENT_BASIC_TYPEDEF(Segment)                                 \
    using SegmentModelDerived = typename traits<Segment>::SegmentModelDerived; \
    using SegmentDataDerived = typename traits<Segment>::SegmentDataDerived;   \
    using NodeDataVector = typename traits<Segment>::NodeDataVector;

#define GALILEO_SEGMENT_CONSTANTS(Segment) \
    static constexpr int NumStages = traits<Segment>::NumStages;

#define GALILEO_SEGMENT_MODEL_TYPEDEF(Segment)                                 \
    using ControlParamModel_t = typename traits<Segment>::ControlParamModel_t; \
    using StageCoefficients_t = typename traits<Segment>::StageCoefficients_t; \
    using Quadrature_t = typename traits<Segment>::Quadrature_t;               \
    using Timings_t = typename traits<Segment>::Timings_t;

#define GALILEO_SEGMENT_DATA_TYPEDEF(Segment)                                \
    using ControlParamData_t = typename traits<Segment>::ControlParamData_t; \
    using S_t = typename traits<Segment>::S_t;                               \
    using Sk_t = typename traits<Segment>::Sk_t;                             \
    using Sw_t = typename traits<Segment>::Sw_t;

namespace galileo
{
    namespace predictive
    {

        template <typename Derived, typename PhaseSpec>
        class SegmentModelBase : internal::CRTP<Derived>
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
            GALILEO_SEGMENT_MODEL_TYPEDEF(SegmentDerived);

            template <typename StateMatrixType, typename ControlMatrixType>
            void calc(SegmentDataDerived &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calc(data, xs.derived(), us.derived());
            }

            template <typename StateMatrixType, typename ControlMatrixType>
            void calcDiff(SegmentDataDerived &data,
                          const Eigen::MatrixBase<StateMatrixType> &xs,
                          const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calcDiff(data, xs.derived(), us.derived());
            }

            template <typename StateMatrixType, typename ControlMatrixType>
            void quasiStatic(SegmentDataDerived &data, const Eigen::MatrixBase<StateMatrixType> &xs,
                             Eigen::MatrixBase<ControlMatrixType> &us,
                             const std::size_t maxiter, const Scalar tol) const
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

            NodeModelDerived node_;

            StageCoefficients_t stage_coefficients_;
            Quadrature_t quadrature_;
            Timings_t timings_;
            Scalar period_;

            ControlParamModel_t *control_parameterization_;
            State_t *state_;
            ActuationModel_t *actuation_;
            ConstraintModelCollection_t *constraints_;
            CostModelCollection_t *costs_;

        }; // class SegmentModelBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_model_base_hpp__
