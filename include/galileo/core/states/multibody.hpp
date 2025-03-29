#ifndef __galileo_core_states_multibody_hpp__
#define __galileo_core_states_multibody_hpp__

#include "galileo/core/states/state-base.hpp"

#include <pinocchio/multibody/model.hpp>

namespace galileo
{

    namespace core
    {

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  int NQ,
                  int NV,
                  int NFb>
        struct StateMultibodyTpl;

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NQ,
                  int _NV,
                  int _NFb>
        struct traits<StateMultibodyTpl<_VarScalar, _NumScalar, _Options, _NQ, _NV, _NFb>>
        {
            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NQ = _NQ;
            static constexpr int NV = _NV;
            static constexpr int NFb = _NFb;

            static constexpr int NX = NQ + NV;
            static constexpr int NU = NV - NFb;
            static constexpr int NDX = NV + NV;

            using VectorNFb_t = Eigen::Matrix<VarScalar, NFb, 1, Options>;
            using VectorNX_t = Eigen::Matrix<VarScalar, NX, 1, Options>;
            using VectorNU_t = Eigen::Matrix<VarScalar, NU, 1, Options>;
            using VectorNDX_t = Eigen::Matrix<VarScalar, NDX, 1, Options>;
            using MatrixNDX_t = Eigen::Matrix<VarScalar, NDX, NDX, Options>;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NQ,
                  int _NV,
                  int _NFb>
        class StateMultibodyTpl : public StateBase<StateMultibodyTpl<_VarScalar, _NumScalar, _Options, _NQ, _NV, _NFb>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using StateDerived = StateMultibodyTpl<_VarScalar, _NumScalar, _Options, _NQ, _NV, _NFb>;
            GALILEO_STATE_BASIC_TYPEDEF(StateDerived);
            GALILEO_STATE_CONSTANTS(StateDerived);
            GALILEO_STATE_TYPEDEF(StateDerived);

            StateMultibodyTpl(
                pinocchio::ModelTpl<VarScalar> *model)
                : model_(model),
                  x0_(VectorNX_t::Zero(NX))
            {
                x0_.head(NQ) = pinocchio::neutral(*model_);

                lb_.head(NFb) =
                    -std::numeric_limits<VarScalar>::infinity() * VectorNFb_t::Ones(NFb);
                ub_.head(NFb) = std::numeric_limits<VarScalar>::infinity() * VectorNFb_t::Ones(NFb);
                lb_.segment(NFb, NQ - NFb) = pinocchio_->lowerPositionLimit.tail(NQ - NFb);
                ub_.segment(NFb, NQ - NFb) = pinocchio_->upperPositionLimit.tail(NQ - NFb);
                lb_.tail(NV) = -pinocchio_->velocityLimit;
                ub_.tail(NV) = pinocchio_->velocityLimit;
            }

            StateMultibodyTpl()
                : x0_(VectorNX_t::Zero(NX)) {}

            ~StateMultibodyTpl() {}

            VectorNX_t zero() const
            {
                return x0_;
            }

            VectorNX_t rand() const
            {
                VectorNX_t xrand = VectorNX_t::Random(NX);
                xrand.head(NQ) = pinocchio::randomConfiguration(*model_);
                return xrand;
            }

            template <typename StateVector1, typename StateVector2, typename StateTangentVector>
            void diff(const Eigen::MatrixBase<StateVector1> &x0,
                      const Eigen::MatrixBase<StateVector2> &x1,
                      Eigen::MatrixBase<StateTangentVector> &dxout) const
            {
                pinocchio::difference(*model_, x0.head(NQ), x1.head(NQ),
                                      dxout.head(NV));
                dxout.tail(NV) = x1.tail(NV) - x0.tail(NV);
            }

            template <typename StateVector1, typename StateTangentVector, typename StateVector2>
            void integrate(const Eigen::MatrixBase<StateVector1> &x,
                           const Eigen::MatrixBase<StateTangentVector> &dx,
                           Eigen::MatrixBase<StateVector2> &xout) const
            {
                pinocchio::integrate(*model_, x.head(NQ), dx.head(NV),
                                     xout.head(NQ));
                xout.tail(NV) = x.tail(NV) + dx.tail(NV);
            }

            template <typename StateVector1, typename StateVector2, typename JMatrix1, typename JMatrix2>
            void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                       const Eigen::MatrixBase<StateVector2> &x1,
                       Eigen::MatrixBase<JMatrix1> &Jfirst, Eigen::MatrixBase<JMatrix2> &Jsecond,
                       const Jcomponent firstsecond = both) const
            {
                if (firstsecond == first)
                {
                    pinocchio::dDifference(*model_, x0.head(NQ), x1.head(NQ),
                                           Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0);
                    Jfirst.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)-1;
                }
                else if (firstsecond == second)
                {
                    pinocchio::dDifference(*model_, x0.head(NQ), x1.head(NQ),
                                           Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1);
                    Jsecond.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)1;
                }
                else
                { // computing both
                    pinocchio::dDifference(*model_, x0.head(NQ), x1.head(NQ),
                                           Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0);
                    pinocchio::dDifference(*model_, x0.head(NQ), x1.head(NQ),
                                           Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1);
                    Jfirst.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)-1;
                    Jsecond.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)1;
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
                        pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                              Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0,
                                              pinocchio::SETTO);
                        Jfirst.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)1;
                        break;
                    case addto:
                        pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                              Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0,
                                              pinocchio::ADDTO);
                        Jfirst.bottomRightCorner(NV, NV).diagonal().array() += (VarScalar)1;
                        break;
                    case rmfrom:
                        pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                              Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0,
                                              pinocchio::RMTO);
                        Jfirst.bottomRightCorner(NV, NV).diagonal().array() -= (VarScalar)1;
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
                        pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                              Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1,
                                              pinocchio::SETTO);
                        Jsecond.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)1;
                        break;
                    case addto:
                        pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                              Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1,
                                              pinocchio::ADDTO);
                        Jsecond.bottomRightCorner(NV, NV).diagonal().array() += (VarScalar)1;
                        break;
                    case rmfrom:
                        pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                              Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1,
                                              pinocchio::RMTO);
                        Jsecond.bottomRightCorner(NV, NV).diagonal().array() -= (VarScalar)1;
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
                    pinocchio::dIntegrateTransport(*model_, x.head(NQ),
                                                   dx.head(NV), Jin.topRows(NV),
                                                   pinocchio::ARG0);
                    break;
                case second:
                    pinocchio::dIntegrateTransport(*model_, x.head(NQ),
                                                   dx.head(NV), Jin.topRows(NV),
                                                   pinocchio::ARG1);
                    break;
                default:
                    break;
                }
            }

            const pinocchio::ModelTpl<VarScalar> *get_model() const
            {
                return model_;
            }

        protected:
            pinocchio::ModelTpl<VarScalar> *model_;
            VectorNX_t x0_;

        }; // class StateMultibodyTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_states_multibody_hpp__