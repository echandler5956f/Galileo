#ifndef __galileo_core_states_multibody_hpp__
#define __galileo_core_states_multibody_hpp__

#include "galileo/core/states/state-base.hpp"

#include <pinocchio/multibody/model.hpp>

namespace galileo
{

    namespace core
    {

        template <typename PhaseSpec, int _NFbq>
        class StateMultibodyTpl : public StateBase<StateMultibodyTpl<PhaseSpec, _NFbq>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
            static constexpr int NFbq = _NFbq;
            using VectorNFbq_t = Eigen::Matrix<typename PS::NumScalar, NFbq, 1>;

            StateMultibodyTpl(
                pinocchio::ModelTpl<typename PS::VarScalar> *model)
                : model_(model),
                  x0_(typename PS::VectorNx_t::Zero(PS::NX))
            {
                x0_.head(PS::NQ) = pinocchio::neutral(*model_);

                lb_.head(NFbq) =
                    -std::numeric_limits<typename PS::NumScalar>::infinity() * VectorNFbq_t::Ones(NFbq);
                ub_.head(NFbq) = std::numeric_limits<typename PS::NumScalar>::infinity() * VectorNFbq_t::Ones(NFbq);
                lb_.segment(NFbq, PS::NQ - NFbq) = pinocchio_->lowerPositionLimit.tail(PS::NQ - NFbq);
                ub_.segment(NFbq, PS::NQ - NFbq) = pinocchio_->upperPositionLimit.tail(PS::NQ - NFbq);
                lb_.tail(PS::NV) = -pinocchio_->velocityLimit;
                ub_.tail(PS::NV) = pinocchio_->velocityLimit;
            }

            StateMultibodyTpl()
                : x0_(typename PS::VectorNx_t::Zero(PS::NX)) {}

            ~StateMultibodyTpl() {}

            typename PS::VectorNx_t zero() const
            {
                return x0_;
            }

            typename PS::VectorNx_t rand() const
            {
                typename PS::VectorNx_t xrand = typename PS::VectorNx_t::Random(PS::NX);
                xrand.head(PS::NQ) = pinocchio::randomConfiguration(*model_);
                return xrand;
            }

            template <typename StateVector1, typename StateVector2, typename StateTangentVector>
            void diff(const Eigen::MatrixBase<StateVector1> &x0,
                      const Eigen::MatrixBase<StateVector2> &x1,
                      Eigen::MatrixBase<StateTangentVector> &dxout) const
            {
                pinocchio::difference(*model_, x0.head(PS::NQ), x1.head(PS::NQ),
                                      dxout.head(PS::NV));
                dxout.tail(PS::NV) = x1.tail(PS::NV) - x0.tail(PS::NV);
            }

            template <typename StateVector1, typename StateTangentVector, typename StateVector2>
            void integrate(const Eigen::MatrixBase<StateVector1> &x,
                           const Eigen::MatrixBase<StateTangentVector> &dx,
                           Eigen::MatrixBase<StateVector2> &xout) const
            {
                pinocchio::integrate(*model_, x.head(PS::NQ), dx.head(PS::NV),
                                     xout.head(PS::NQ));
                xout.tail(PS::NV) = x.tail(PS::NV) + dx.tail(PS::NV);
            }

            template <typename StateVector1, typename StateVector2, typename JMatrix1, typename JMatrix2>
            void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                       const Eigen::MatrixBase<StateVector2> &x1,
                       Eigen::MatrixBase<JMatrix1> &Jfirst, Eigen::MatrixBase<JMatrix2> &Jsecond,
                       const Jcomponent firstsecond = both) const
            {
                if (firstsecond == first)
                {
                    pinocchio::dDifference(*model_, x0.head(PS::NQ), x1.head(PS::NQ),
                                           Jfirst.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG0);
                    Jfirst.bottomRightCorner(PS::NV, PS::NV).diagonal().array() = (typename PS::VarScalar) - 1;
                }
                else if (firstsecond == second)
                {
                    pinocchio::dDifference(*model_, x0.head(PS::NQ), x1.head(PS::NQ),
                                           Jsecond.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG1);
                    Jsecond.bottomRightCorner(PS::NV, PS::NV).diagonal().array() = typename PS::VarScalar(1);
                }
                else
                { // computing both
                    pinocchio::dDifference(*model_, x0.head(PS::NQ), x1.head(PS::NQ),
                                           Jfirst.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG0);
                    pinocchio::dDifference(*model_, x0.head(PS::NQ), x1.head(PS::NQ),
                                           Jsecond.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG1);
                    Jfirst.bottomRightCorner(PS::NV, PS::NV).diagonal().array() = (typename PS::VarScalar) - 1;
                    Jsecond.bottomRightCorner(PS::NV, PS::NV).diagonal().array() = (typename PS::VarScalar)1;
                }
            }

            template <typename StateVector, typename StateTangentVector, typename JMatrix1, typename JMatrix2>
            void Jintegrate(const Eigen::MatrixBase<StateVector> &x,
                            const Eigen::MatrixBase<StateTangentVector> &dx,
                            Eigen::MatrixBase<JMatrix1> &Jfirst,
                            Eigen::MatrixBase<JMatrix2> &Jsecond,
                            const Jcomponent firstsecond = both,
                            const AssignmentOp op = setto) const
            {
                if (firstsecond == first || firstsecond == both)
                {
                    switch (op)
                    {
                    case setto:
                        pinocchio::dIntegrate(*model_, x.head(PS::NQ), dx.head(PS::NV),
                                              Jfirst.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG0,
                                              pinocchio::SETTO);
                        Jfirst.bottomRightCorner(PS::NV, PS::NV).diagonal().array() = typename PS::VarScalar(1);
                        break;
                    case addto:
                        pinocchio::dIntegrate(*model_, x.head(PS::NQ), dx.head(PS::NV),
                                              Jfirst.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG0,
                                              pinocchio::ADDTO);
                        Jfirst.bottomRightCorner(PS::NV, PS::NV).diagonal().array() += typename PS::VarScalar(1);
                        break;
                    case rmfrom:
                        pinocchio::dIntegrate(*model_, x.head(PS::NQ), dx.head(PS::NV),
                                              Jfirst.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG0,
                                              pinocchio::RMTO);
                        Jfirst.bottomRightCorner(PS::NV, PS::NV).diagonal().array() -= typename PS::VarScalar(1);
                        break;
                    default:
                        break;
                    }
                }
                if (firstsecond == second || firstsecond == both)
                {
                    switch (op)
                    {
                    case setto:
                        pinocchio::dIntegrate(*model_, x.head(PS::NQ), dx.head(PS::NV),
                                              Jsecond.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG1,
                                              pinocchio::SETTO);
                        Jsecond.bottomRightCorner(PS::NV, PS::NV).diagonal().array() = typename PS::VarScalar(1);
                        break;
                    case addto:
                        pinocchio::dIntegrate(*model_, x.head(PS::NQ), dx.head(PS::NV),
                                              Jsecond.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG1,
                                              pinocchio::ADDTO);
                        Jsecond.bottomRightCorner(PS::NV, PS::NV).diagonal().array() += typename PS::VarScalar(1);
                        break;
                    case rmfrom:
                        pinocchio::dIntegrate(*model_, x.head(PS::NQ), dx.head(PS::NV),
                                              Jsecond.topLeftCorner(PS::NV, PS::NV), pinocchio::ARG1,
                                              pinocchio::RMTO);
                        Jsecond.bottomRightCorner(PS::NV, PS::NV).diagonal().array() -= typename PS::VarScalar(1);
                        break;
                    default:
                        break;
                    }
                }
            }

            template <typename StateVector, typename StateTangentVector, typename JMatrix>
            void JintegrateTransport(const Eigen::MatrixBase<StateVector> &x,
                                     const Eigen::MatrixBase<StateTangentVector> &dx,
                                     Eigen::MatrixBase<JMatrix> &Jin,
                                     const Jcomponent firstsecond) const
            {
                switch (firstsecond)
                {
                case first:
                    pinocchio::dIntegrateTransport(*model_, x.head(PS::NQ),
                                                   dx.head(PS::NV), Jin.topRows(PS::NV),
                                                   pinocchio::ARG0);
                    break;
                case second:
                    pinocchio::dIntegrateTransport(*model_, x.head(PS::NQ),
                                                   dx.head(PS::NV), Jin.topRows(PS::NV),
                                                   pinocchio::ARG1);
                    break;
                default:
                    break;
                }
            }

            const pinocchio::ModelTpl<typename PS::VarScalar> *get_model() const
            {
                return model_;
            }

        protected:
            pinocchio::ModelTpl<typename PS::VarScalar> *model_;
            typename PS::VectorNx_t x0_;

        }; // class StateMultibodyTpl

        // A wrapper that collapses <PS, _NFbq> into a single template <PS>.
        template <int _NFbq>
        struct StateMetaMultibody
        {
            template <typename PhaseSpec>
            struct Implementation
            {
                using State_t = StateMultibodyTpl<PhaseSpec, _NFbq>;
            };
        };

    } // namespace core

} // namespace galileo

#endif // __galileo_core_states_multibody_hpp__