#ifndef __galileo_predictive_phases_phase_model_base_hpp__
#define __galileo_predictive_phases_phase_model_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/segments/segment-base.hpp"

#define GALILEO_PHASE_BASIC_TYPEDEF(Phase)                                 \
    using PhaseModelDerived = typename traits<Phase>::PhaseModelDerived;   \
    using PhaseDataDerived = typename traits<Phase>::PhaseDataDerived;     \
    using SegmentModelVector = typename traits<Phase>::SegmentModelVector; \
    using SegmentDataVector = typename traits<Phase>::SegmentDataVector;

#define GALILEO_PHASE_CONSTANTS(Phase) \
    static constexpr int NumSegments = traits<Phase>::NumSegments;

#define GALILEO_PHASE_MODEL_TYPEDEF(Phase)

#define GALILEO_PHASE_DATA_TYPEDEF(Phase)

namespace galileo
{
    namespace predictive
    {

        template <typename Derived>
        class PhaseModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using NodeDerived = typename traits<Derived>::NodeDerived;
            using SegmentDerived = typename traits<Derived>::SegmentDerived;
            using PhaseDerived = typename traits<Derived>::PhaseDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_BASIC_TYPEDEF(SegmentDerived);
            GALILEO_PHASE_BASIC_TYPEDEF(PhaseDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_CONSTANTS(SegmentDerived);
            GALILEO_PHASE_CONSTANTS(PhaseDerived);

            GALILEO_NODE_MODEL_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_MODEL_TYPEDEF(SegmentDerived);
            GALILEO_PHASE_MODEL_TYPEDEF(PhaseDerived);

            template <typename StateMatrixType, typename ControlMatrixType>
            void calc(PhaseDataDerived &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calc(data, xs.derived(), us.derived());
            }

            template <typename StateMatrixType, typename ControlMatrixType>
            void calcDiff(PhaseDataDerived &data,
                          const Eigen::MatrixBase<StateMatrixType> &xs,
                          const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
                derived().calcDiff(data, xs.derived(), us.derived());
            }

            template <typename StateMatrixType, typename ControlMatrixType>
            void quasiStatic(PhaseDataDerived &data, const Eigen::MatrixBase<StateMatrixType> &xs,
                             Eigen::MatrixBase<ControlMatrixType> &us,
                             const std::size_t maxiter, const Scalar tol) const
            {
                derived().quasiStatic(data, xs.derived(), us.derived(), maxiter, tol);
            }

        protected:
            inline PhaseModelBase()
            {
            }

            inline PhaseModelBase(const PhaseModelBase &clone)
            {
                *this = clone;
            }

            inline PhaseModelBase &operator=(const PhaseModelBase &clone)
            {
                return *this;
            }

            SegmentModelVector segments_;

            Scalar phase_period_;

            ControlParamModel_t control_parameterization_;
            State_t state_;
            ActuationModel_t actuation_;
            ConstraintModelCollection_t constraints_;
            CostModelCollection_t costs_;

        }; // class PhaseModelBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_model_base_hpp__
