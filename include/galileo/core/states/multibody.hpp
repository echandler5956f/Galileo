#ifndef __galileo_core_states_multibody_hpp__
#define __galileo_core_states_multibody_hpp__

#include <pinocchio/multibody/model.hpp>

#include "galileo/core/states/state-base.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  int NX,
                  int NU,
                  int NDX>
        struct StateMultibodyTpl;

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NX,
                  int _NU,
                  int _NDX>
        struct traits<StateMultibodyTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>>
        {

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NX = _NX;
            static constexpr int NU = _NU;
            static constexpr int NDX = _NDX;

            using VectorNX_t = Eigen::Matrix<VarScalar, NX, 1, Options>;
            using VectorNU_t = Eigen::Matrix<VarScalar, NU, 1, Options>;
            using VectorNDX_t = Eigen::Matrix<VarScalar, NDX, 1, Options>;
            using MatrixNDX_t = Eigen::Matrix<VarScalar, NDX, NDX, Options>;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NX,
                  int _NU,
                  int _NDX>
        class StateMultibodyTpl : public StateBase<StateMultibodyTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using StateDerived = StateMultibodyTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>;
            GALILEO_STATE_BASIC_TYPEDEF(StateDerived);
            GALILEO_STATE_CONSTANTS(StateDerived);
            GALILEO_STATE_TYPEDEF(StateDerived);

            StateMultibodyTpl(
                pinocchio::ModelTpl<VarScalar> *model)
                : model_(model),
                  x0_(VectorNX_t::Zero(model->nq + model->nv)),
                  nq_(model->nq),
                  nv_(model->nv)
            {
                x0_.head(nq_) = pinocchio::neutral(*model_);

                const std::size_t nq0 = model->joints[1].nq();

                lb_.head(nq0) =
                    -std::numeric_limits<VarScalar>::infinity() * VectorNX_t::Ones(nq0);
                ub_.head(nq0) = std::numeric_limits<VarScalar>::infinity() * VectorNX_t::Ones(nq0);
                lb_.segment(nq0, nq_ - nq0) = pinocchio_->lowerPositionLimit.tail(nq_ - nq0);
                ub_.segment(nq0, nq_ - nq0) = pinocchio_->upperPositionLimit.tail(nq_ - nq0);
                lb_.tail(nv_) = -pinocchio_->velocityLimit;
                ub_.tail(nv_) = pinocchio_->velocityLimit;
            }

            StateMultibodyTpl()
                : x0_(VectorNX_t::Zero(0)) {}

            ~StateMultibodyTpl() {}

            VectorNX_t zero() const
            {
                return x0_;
            }

            VectorNX_t rand() const
            {
                VectorNX_t xrand = VectorNX_t::Random(nx_);
                xrand.head(nq_) = pinocchio::randomConfiguration(*model_);
                return xrand;
            }

            void diff(const Eigen::Ref<const VectorNX_t> &x0,
                      const Eigen::Ref<const VectorNX_t> &x1,
                      Eigen::Ref<VectorNX_t> dxout) const
            {
                pinocchio::difference(*model_, x0.head(nq_), x1.head(nq_),
                                      dxout.head(nv_));
                dxout.tail(nv_) = x1.tail(nv_) - x0.tail(nv_);
            }

            void integrate(const Eigen::Ref<const VectorNX_t> &x,
                           const Eigen::Ref<const VectorNX_t> &dx,
                           Eigen::Ref<VectorNX_t> xout) const
            {
                pinocchio::integrate(*model_, x.head(nq_), dx.head(nv_),
                                     xout.head(nq_));
                xout.tail(nv_) = x.tail(nv_) + dx.tail(nv_);
            }

            void Jdiff(const Eigen::Ref<const VectorNX_t> &x0,
                       const Eigen::Ref<const VectorNX_t> &x1,
                       Eigen::Ref<MatrixNDX_t> Jfirst,
                       Eigen::Ref<MatrixNDX_t> Jsecond,
                       const Jcomponent firstsecond) const
            {
                if (firstsecond == first)
                {
                    pinocchio::dDifference(*model_, x0.head(nq_), x1.head(nq_),
                                           Jfirst.topLeftCorner(nv_, nv_), pinocchio::ARG0);
                    Jfirst.bottomRightCorner(nv_, nv_).diagonal().array() = (VarScalar)-1;
                }
                else if (firstsecond == second)
                {
                    pinocchio::dDifference(*model_, x0.head(nq_), x1.head(nq_),
                                           Jsecond.topLeftCorner(nv_, nv_), pinocchio::ARG1);
                    Jsecond.bottomRightCorner(nv_, nv_).diagonal().array() = (VarScalar)1;
                }
                else
                { // computing both
                    pinocchio::dDifference(*model_, x0.head(nq_), x1.head(nq_),
                                           Jfirst.topLeftCorner(nv_, nv_), pinocchio::ARG0);
                    pinocchio::dDifference(*model_, x0.head(nq_), x1.head(nq_),
                                           Jsecond.topLeftCorner(nv_, nv_), pinocchio::ARG1);
                    Jfirst.bottomRightCorner(nv_, nv_).diagonal().array() = (VarScalar)-1;
                    Jsecond.bottomRightCorner(nv_, nv_).diagonal().array() = (VarScalar)1;
                }
            }

            void Jintegrate(const Eigen::Ref<const VectorNX_t> &x,
                            const Eigen::Ref<const VectorNX_t> &dx,
                            Eigen::Ref<MatrixNDX_t> Jfirst,
                            Eigen::Ref<MatrixNDX_t> Jsecond,
                            const Jcomponent firstsecond,
                            const AssignmentOp op) const
            {
                if (firstsecond == first || firstsecond == both)
                {
                    switch (op)
                    {
                    case setto:
                        pinocchio::dIntegrate(*model_, x.head(nq_), dx.head(nv_),
                                              Jfirst.topLeftCorner(nv_, nv_), pinocchio::ARG0,
                                              pinocchio::SETTO);
                        Jfirst.bottomRightCorner(nv_, nv_).diagonal().array() = (VarScalar)1;
                        break;
                    case addto:
                        pinocchio::dIntegrate(*model_, x.head(nq_), dx.head(nv_),
                                              Jfirst.topLeftCorner(nv_, nv_), pinocchio::ARG0,
                                              pinocchio::ADDTO);
                        Jfirst.bottomRightCorner(nv_, nv_).diagonal().array() += (VarScalar)1;
                        break;
                    case rmfrom:
                        pinocchio::dIntegrate(*model_, x.head(nq_), dx.head(nv_),
                                              Jfirst.topLeftCorner(nv_, nv_), pinocchio::ARG0,
                                              pinocchio::RMTO);
                        Jfirst.bottomRightCorner(nv_, nv_).diagonal().array() -= (VarScalar)1;
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
                        pinocchio::dIntegrate(*model_, x.head(nq_), dx.head(nv_),
                                              Jsecond.topLeftCorner(nv_, nv_), pinocchio::ARG1,
                                              pinocchio::SETTO);
                        Jsecond.bottomRightCorner(nv_, nv_).diagonal().array() = (VarScalar)1;
                        break;
                    case addto:
                        pinocchio::dIntegrate(*model_, x.head(nq_), dx.head(nv_),
                                              Jsecond.topLeftCorner(nv_, nv_), pinocchio::ARG1,
                                              pinocchio::ADDTO);
                        Jsecond.bottomRightCorner(nv_, nv_).diagonal().array() += (VarScalar)1;
                        break;
                    case rmfrom:
                        pinocchio::dIntegrate(*model_, x.head(nq_), dx.head(nv_),
                                              Jsecond.topLeftCorner(nv_, nv_), pinocchio::ARG1,
                                              pinocchio::RMTO);
                        Jsecond.bottomRightCorner(nv_, nv_).diagonal().array() -= (VarScalar)1;
                        break;
                    default:
                        break;
                    }
                }
            }

            void JintegrateTransport(const Eigen::Ref<const VectorNX_t> &x,
                                     const Eigen::Ref<const VectorNX_t> &dx,
                                     Eigen::Ref<MatrixNDX_t> Jin,
                                     const Jcomponent firstsecond) const
            {
                switch (firstsecond)
                {
                case first:
                    pinocchio::dIntegrateTransport(*model_, x.head(nq_),
                                                   dx.head(nv_), Jin.topRows(nv_),
                                                   pinocchio::ARG0);
                    break;
                case second:
                    pinocchio::dIntegrateTransport(*model_, x.head(nq_),
                                                   dx.head(nv_), Jin.topRows(nv_),
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

            std::size_t nq_;
            std::size_t nv_;

        }; // class StateMultibodyTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_states_multibody_hpp__