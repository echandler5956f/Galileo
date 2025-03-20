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

            NumScalar period() const
            {
                return derived().period();
            }

            template <typename StateVectorType, typename ControlVectorType>
            void segmentCalc(const std::size_t &segment_index,
                             const Eigen::MatrixBase<StateVectorType> &xs,
                             const Eigen::MatrixBase<ControlVectorType> &us) const
            {
                derived().segmentCalc(segment_index, xs.derived(), us.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void segmentCalcDiff(const std::size_t &segment_index,
                                 const Eigen::MatrixBase<StateVectorType> &xs,
                                 const Eigen::MatrixBase<ControlVectorType> &us) const
            {
                derived().segmentCalcDiff(segment_index, xs.derived(), us.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void segmentQuasiStatic(const std::size_t &segment_index,
                                    const Eigen::MatrixBase<StateVectorType> &xs,
                                    Eigen::MatrixBase<ControlVectorType> &us,
                                    const std::size_t maxiter, const Scalar tol) const
            {
                derived().segmentQuasiStatic(segment_index, xs.derived(), us.derived(), maxiter, tol);
            }

            State_t::VectorNX_t stateZero() const
            {
                return derived().stateZero();
            }

            State_t::VectorNX_t stateRand() const
            {
                return derived().stateRand();
            }

            template <typename StateVectorType1, typename StateVectorType2, typename StateTangentVectorType>
            void stateDiff(const Eigen::MatrixBase<StateVectorType1> &x0, const Eigen::MatrixBase<StateVectorType2> &x1, Eigen::MatrixBase<StateTangentVectorType> &dxout) const
            {
                derived().stateDiff(x0.derived(), x1.derived(), dxout.derived());
            }

            template <typename StateVectorType, typename StateTangentVectorType, typename StateVectorType2>
            void stateIntegrate(const Eigen::MatrixBase<StateVectorType> &x, const Eigen::MatrixBase<StateTangentVectorType> &dx, Eigen::MatrixBase<StateVectorType2> &xout) const
            {
                derived().stateIntegrate(x.derived(), dx.derived(), xout.derived());
            }

            template <typename StateVectorType1, typename StateTangentVectorType, typename StateVectorType2>
            void stateJdiff(const Eigen::MatrixBase<StateVectorType1> &x0, const Eigen::MatrixBase<StateVectorType2> &x1, Eigen::MatrixBase<JMatrix1> &Jfirst, Eigen::MatrixBase<JMatrix2> &Jsecond, const Jcomponent firstsecond = Jcomponent::both) const
            {
                derived().stateJdiff(x0.derived(), x1.derived(), Jfirst.derived(), Jsecond.derived(), firstsecond);
            }

            template <typename StateVectorType, typename StateTangentVectorType, typename StateVectorType2>
            void stateJintegrate(const Eigen::MatrixBase<StateVectorType> &x, const Eigen::MatrixBase<StateTangentVectorType> &dx, Eigen::MatrixBase<JMatrix1> &Jfirst, Eigen::MatrixBase<JMatrix2> &Jsecond, const Jcomponent firstsecond = Jcomponent::both, const AssignmentOp op = AssignmentOp::setto) const
            {
                derived().stateJintegrate(x.derived(), dx.derived(), Jfirst.derived(), Jsecond.derived(), firstsecond, op);
            }

            template <typename StateVectorType, typename StateTangentVectorType, typename JMatrix>
            void stateJintegrateTransport(const Eigen::MatrixBase<StateVectorType> &x, const Eigen::MatrixBase<StateTangentVectorType> &dx, Eigen::MatrixBase<JMatrix> &Jin, const Jcomponent firstsecond) const
            {
                derived().stateJintegrateTransport(x.derived(), dx.derived(), Jin.derived(), firstsecond);
            }

            template <typename StateVectorType1, typename StateVectorType2, typename StateTangentVectorType>
            State_t::VectorNDX_t stateDiffDx(const Eigen::MatrixBase<StateVectorType1> &x0, const Eigen::MatrixBase<StateVectorType2> &x1) const
            {
                return derived().stateDiffDx(x0.derived(), x1.derived());
            }

            template <typename StateVectorType, typename StateTangentVectorType>
            State_t::VectorNX_t stateIntegrateX(const Eigen::MatrixBase<StateVectorType> &x, const Eigen::MatrixBase<StateTangentVectorType> &dx) const
            {
                return derived().stateIntegrateX(x.derived(), dx.derived());
            }

            template <typename StateVectorType1, typename StateVectorType2>
            std::vector<State_t::MatrixNDX_t> stateJdiffJs(const Eigen::MatrixBase<StateVectorType1> &x0, const Eigen::MatrixBase<StateVectorType2> &x1, const Jcomponent firstsecond = Jcomponent::both) const
            {
                return derived().stateJdiffJs(x0.derived(), x1.derived(), firstsecond);
            }

            template <typename StateVectorType, typename StateTangentVectorType>
            std::vector<State_t::MatrixNDX_t> stateJintegrateJs(const Eigen::MatrixBase<StateVectorType> &x, const Eigen::MatrixBase<StateTangentVectorType> &dx, const Jcomponent firstsecond = Jcomponent::both) const
            {
                return derived().stateJintegrateJs(x.derived(), dx.derived(), firstsecond);
            }

            int nx() const
            {
                return derived().nx();
            }

            int nu() const
            {
                return derived().nu();
            }

            int ndx() const
            {
                return derived().ndx();
            }

            int nh() const
            {
                return derived().nh();
            }

            int ng() const
            {
                return derived().ng();
            }

            int nc() const
            {
                return derived().nc();
            }

            const typename ConstraintModelCollection_t::H_Equality_t &h_eq() const
            {
                return derived().h_eq();
            }

            const typename ConstraintModelCollection_t::G_Bound_t &g_lb() const
            {
                return derived().g_lb();
            }

            const typename ConstraintModelCollection_t::G_Bound_t &g_ub() const
            {
                return derived().g_ub();
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
