#ifndef __galileo_core_states_multibody_hpp__
#define __galileo_core_states_multibody_hpp__

#include "galileo/core/states/state-base.hpp"

#include <pinocchio/multibody/model.hpp>

namespace galileo
{

    namespace core
    {

        template <typename BasicSpec>
        class StateMultibodyTpl : public StateBase<StateMultibodyTpl<BasicSpec>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using BS = BasicSpec;

            StateMultibodyTpl(
                pinocchio::ModelTpl<typename BS::VarScalar> *model)
                : model_(model),
                  x0_(typename BS::VectorNx_t::Zero(BS::NX))
            {
                x0_.head(BS::NQ) = pinocchio::neutral(*model_);

                lb_.head(BS::NQb) =
                    -std::numeric_limits<typename BS::NumScalar>::infinity() * VectorNQb_t::Ones(BS::NQb);
                ub_.head(BS::NQb) = std::numeric_limits<typename BS::NumScalar>::infinity() * VectorNQb_t::Ones(NQb);
                lb_.segment(BS::NQb, BS::NQ - BS::NQb) = pinocchio_->lowerPositionLimit.tail(BS::NQ - BS::NQb);
                ub_.segment(BS::NQb, BS::NQ - BS::NQb) = pinocchio_->upperPositionLimit.tail(BS::NQ - BS::NQb);
                lb_.tail(BS::NV) = -pinocchio_->velocityLimit;
                ub_.tail(BS::NV) = pinocchio_->velocityLimit;
            }

            StateMultibodyTpl()
                : x0_(typename BS::VectorNx_t::Zero(BS::NX)) {}

            ~StateMultibodyTpl() {}

            typename BS::VectorNx_t zero() const
            {
                return x0_;
            }

            typename BS::VectorNx_t rand() const
            {
                typename BS::VectorNx_t xrand = typename BS::VectorNx_t::Random(BS::NX);
                xrand.head(BS::NQ) = pinocchio::randomConfiguration(*model_);
                return xrand;
            }

            template <typename StateVector1, typename StateVector2, typename StateTangentVector>
            void diff(const Eigen::MatrixBase<StateVector1> &x0,
                      const Eigen::MatrixBase<StateVector2> &x1,
                      Eigen::MatrixBase<StateTangentVector> &dxout) const
            {
                pinocchio::difference(*model_, x0.head(BS::NQ), x1.head(BS::NQ),
                                      dxout.head(BS::NV));
                dxout.tail(BS::NV) = x1.tail(BS::NV) - x0.tail(BS::NV);
            }

            template <typename StateVector1, typename StateTangentVector, typename StateVector2>
            void integrate(const Eigen::MatrixBase<StateVector1> &x,
                           const Eigen::MatrixBase<StateTangentVector> &dx,
                           Eigen::MatrixBase<StateVector2> &xout) const
            {
                pinocchio::integrate(*model_, x.head(BS::NQ), dx.head(BS::NV),
                                     xout.head(BS::NQ));
                xout.tail(BS::NV) = x.tail(BS::NV) + dx.tail(BS::NV);
            }

            template <typename StateVector1, typename StateVector2, typename JMatrix1, typename JMatrix2>
            void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                       const Eigen::MatrixBase<StateVector2> &x1,
                       Eigen::MatrixBase<JMatrix1> &Jfirst, Eigen::MatrixBase<JMatrix2> &Jsecond,
                       const Jcomponent firstsecond = both) const
            {
                if (firstsecond == first)
                {
                    pinocchio::dDifference(*model_, x0.head(BS::NQ), x1.head(BS::NQ),
                                           Jfirst.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG0);
                    Jfirst.bottomRightCorner(BS::NV, BS::NV).diagonal().array() = (typename BS::VarScalar) - 1;
                }
                else if (firstsecond == second)
                {
                    pinocchio::dDifference(*model_, x0.head(BS::NQ), x1.head(BS::NQ),
                                           Jsecond.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG1);
                    Jsecond.bottomRightCorner(BS::NV, BS::NV).diagonal().array() = typename BS::VarScalar(1);
                }
                else
                { // computing both
                    pinocchio::dDifference(*model_, x0.head(BS::NQ), x1.head(BS::NQ),
                                           Jfirst.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG0);
                    pinocchio::dDifference(*model_, x0.head(BS::NQ), x1.head(BS::NQ),
                                           Jsecond.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG1);
                    Jfirst.bottomRightCorner(BS::NV, BS::NV).diagonal().array() = (typename BS::VarScalar) - 1;
                    Jsecond.bottomRightCorner(BS::NV, BS::NV).diagonal().array() = (typename BS::VarScalar)1;
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
                        pinocchio::dIntegrate(*model_, x.head(BS::NQ), dx.head(BS::NV),
                                              Jfirst.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG0,
                                              pinocchio::SETTO);
                        Jfirst.bottomRightCorner(BS::NV, BS::NV).diagonal().array() = typename BS::VarScalar(1);
                        break;
                    case addto:
                        pinocchio::dIntegrate(*model_, x.head(BS::NQ), dx.head(BS::NV),
                                              Jfirst.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG0,
                                              pinocchio::ADDTO);
                        Jfirst.bottomRightCorner(BS::NV, BS::NV).diagonal().array() += typename BS::VarScalar(1);
                        break;
                    case rmfrom:
                        pinocchio::dIntegrate(*model_, x.head(BS::NQ), dx.head(BS::NV),
                                              Jfirst.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG0,
                                              pinocchio::RMTO);
                        Jfirst.bottomRightCorner(BS::NV, BS::NV).diagonal().array() -= typename BS::VarScalar(1);
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
                        pinocchio::dIntegrate(*model_, x.head(BS::NQ), dx.head(BS::NV),
                                              Jsecond.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG1,
                                              pinocchio::SETTO);
                        Jsecond.bottomRightCorner(BS::NV, BS::NV).diagonal().array() = typename BS::VarScalar(1);
                        break;
                    case addto:
                        pinocchio::dIntegrate(*model_, x.head(BS::NQ), dx.head(BS::NV),
                                              Jsecond.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG1,
                                              pinocchio::ADDTO);
                        Jsecond.bottomRightCorner(BS::NV, BS::NV).diagonal().array() += typename BS::VarScalar(1);
                        break;
                    case rmfrom:
                        pinocchio::dIntegrate(*model_, x.head(BS::NQ), dx.head(BS::NV),
                                              Jsecond.topLeftCorner(BS::NV, BS::NV), pinocchio::ARG1,
                                              pinocchio::RMTO);
                        Jsecond.bottomRightCorner(BS::NV, BS::NV).diagonal().array() -= typename BS::VarScalar(1);
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
                    pinocchio::dIntegrateTransport(*model_, x.head(BS::NQ),
                                                   dx.head(BS::NV), Jin.topRows(BS::NV),
                                                   pinocchio::ARG0);
                    break;
                case second:
                    pinocchio::dIntegrateTransport(*model_, x.head(BS::NQ),
                                                   dx.head(BS::NV), Jin.topRows(BS::NV),
                                                   pinocchio::ARG1);
                    break;
                default:
                    break;
                }
            }

            const pinocchio::ModelTpl<typename BS::VarScalar> *get_model() const
            {
                return model_;
            }

        protected:
            pinocchio::ModelTpl<typename BS::VarScalar> *model_;
            typename BS::VectorNx_t x0_;

        }; // class StateMultibodyTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_states_multibody_hpp__