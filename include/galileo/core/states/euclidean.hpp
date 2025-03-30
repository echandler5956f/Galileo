#ifndef __galileo_core_states_euclidean_hpp__
#define __galileo_core_states_euclidean_hpp__

#include "galileo/core/states/state-base.hpp"

namespace galileo
{

    namespace core
    {

        template <typename PhaseSpec>
        class StateEuclideanTpl : public StateBase<StateEuclideanTpl<PhaseSpec>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            typename PS::VectorNx_t zero() const
            {
                return typename PS::VectorNx_t::Zero();
            }

            typename PS::VectorNx_t rand() const
            {
                return typename PS::VectorNx_t::Random();
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
                    Jfirst.diagonal() = typename PS::VectorNdx_t::Constant(PS::NDX, -1.);
                }
                if (firstsecond == second || firstsecond == both)
                {
                    Jsecond.setZero();
                    Jsecond.diagonal() = typename PS::VectorNdx_t::Constant(PS::NDX, 1.);
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
                        Jfirst.diagonal().array() = typename PS::VarScalar(1.);
                        break;
                    case addto:
                        Jfirst.diagonal().array() += typename PS::VarScalar(1.);
                        break;
                    case rmfrom:
                        Jfirst.diagonal().array() -= typename PS::VarScalar(1.);
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
                        Jsecond.diagonal().array() = typename PS::VarScalar(1.);
                        break;
                    case addto:
                        Jsecond.diagonal().array() += typename PS::VarScalar(1.);
                        break;
                    case rmfrom:
                        Jsecond.diagonal().array() -= typename PS::VarScalar(1.);
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