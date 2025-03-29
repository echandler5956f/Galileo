#ifndef __galileo_core_states_euclidean_hpp__
#define __galileo_core_states_euclidean_hpp__

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
        struct StateEuclideanTpl;

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NX,
                  int _NU,
                  int _NDX>
        struct traits<StateEuclideanTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>>
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
        class StateEuclideanTpl : public StateBase<StateEuclideanTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using StateDerived = StateEuclideanTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>;
            GALILEO_STATE_BASIC_TYPEDEF(StateDerived);
            GALILEO_STATE_CONSTANTS(StateDerived);
            GALILEO_STATE_TYPEDEF(StateDerived);

            VectorNX_t zero() const
            {
                return VectorNX_t::Zero();
            }

            VectorNX_t rand() const
            {
                return VectorNX_t::Random();
            }

            template <typename StateVector1, typename StateVector2, typename StateTangentVector>
            void diff(const Eigen::MatrixBase<StateVector1> &x0,
                      const Eigen::MatrixBase<StateVector2> &x1,
                      Eigen::MatrixBase<StateTangentVector> &dxout) const
            {
                dxout = x1.derived() - x0.derived();
            }

            template <typename StateVector1, typename StateTangentVector, typename StateVector2>
            void integrate(const Eigen::MatrixBase<StateVector1> &x,
                           const Eigen::MatrixBase<StateTangentVector> &dx,
                           Eigen::MatrixBase<StateVector2> &xout) const
            {
                xout = x.derived() + dx.derived();
            }

            template <typename StateVector1, typename StateVector2, typename JMatrix1, typename JMatrix2>
            void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                       const Eigen::MatrixBase<StateVector2> &x1,
                       Eigen::MatrixBase<JMatrix1> &Jfirst, Eigen::MatrixBase<JMatrix2> &Jsecond,
                       const Jcomponent firstsecond = both) const
            {
                if (firstsecond == first || firstsecond == both)
                {
                    Jfirst.setZero();
                    Jfirst.diagonal() = VectorNDX_t::Constant(NDX, -1.);
                }
                if (firstsecond == second || firstsecond == both)
                {
                    Jsecond.setZero();
                    Jsecond.diagonal() = VectorNDX_t::Constant(NDX, 1.);
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
                        Jfirst.diagonal().array() = Scalar(1.);
                        break;
                    case addto:
                        Jfirst.diagonal().array() += Scalar(1.);
                        break;
                    case rmfrom:
                        Jfirst.diagonal().array() -= Scalar(1.);
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
                        Jsecond.diagonal().array() = Scalar(1.);
                        break;
                    case addto:
                        Jsecond.diagonal().array() += Scalar(1.);
                        break;
                    case rmfrom:
                        Jsecond.diagonal().array() -= Scalar(1.);
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
                // todo: not implemented
            }

        }; // class StateEuclideanTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_states_euclidean_hpp__