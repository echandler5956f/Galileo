#ifndef __galileo_core_states_euclidean_hpp__
#define __galileo_core_states_euclidean_hpp__

#include "galileo/core/states/state-base.hpp"
#include "galileo/multibody/robot-spec.hpp"

namespace galileo
{

    template <typename RobotSpec>
    class StateEuclideanTpl
        : public StateBase<StateEuclideanTpl<RobotSpec>, RobotSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

        GALILEO_ROBOT_SPEC_SCALARS_TYPEDEF(RS);
        GALILEO_ROBOT_SPEC_EIGEN_TYPES_TYPEDEF(RS);

        using Base = StateBase<StateEuclideanTpl<RS>, RS>;

        // For StateEuclideanTpl specifically, robot spec must be valid and all dimensions must be set apriori.
        StateEuclideanTpl(const RS &rs, const VectorNx_t &lb, const VectorNx_t &ub)
            : Base(rs), lb_(lb), ub_(ub)
        {
            GALILEO_ASSERT(IsValidRobotSpec(rs), "StateEuclideanTpl: Invalid robot spec");
            GALILEO_ASSERT(get_nx() == get_ndx(), "StateEuclideanTpl: Invalid dimensions for Euclidean state (nx != ndx)");
        }

        VectorNx_t zero() const
        {
            return VectorNx_t::Zero(get_nx());
        }

        VectorNx_t rand() const
        {
            return VectorNx_t::Random(get_nx());
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
                Jfirst.diagonal() = VectorNdx_t::Constant(get_ndx(), VarScalar(-1.));
            }
            if (firstsecond == second || firstsecond == both)
            {
                Jsecond.setZero();
                Jsecond.diagonal() = VectorNdx_t::Constant(get_ndx(), VarScalar(1.));
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
                    Jfirst.diagonal().array() = VarScalar(1.);
                    break;
                case addto:
                    Jfirst.diagonal().array() += VarScalar(1.);
                    break;
                case rmfrom:
                    Jfirst.diagonal().array() -= VarScalar(1.);
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
                    Jsecond.diagonal().array() = VarScalar(1.);
                    break;
                case addto:
                    Jsecond.diagonal().array() += VarScalar(1.);
                    break;
                case rmfrom:
                    Jsecond.diagonal().array() -= VarScalar(1.);
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
            // Nothing to do
        }

        const VectorNx_t &get_lb_impl() const
        {
            return lb_;
        }

        const VectorNx_t &get_ub_impl() const
        {
            return ub_;
        }

        template <typename StateVector>
        void set_lb_impl(const Eigen::MatrixBase<StateVector> &lb)
        {
            lb_ = lb.derived();
        }

        template <typename StateVector>
        void set_ub_impl(const Eigen::MatrixBase<StateVector> &ub)
        {
            ub_ = ub.derived();
        }

        using Base::diff_dx;
        using Base::integrate_x;
        using Base::Jdiff_Js;
        using Base::Jintegrate_Js;

        using Base::get_rs;

        using Base::get_nqb;
        using Base::get_nqb_dim;

        using Base::get_nqj;
        using Base::get_nqj_dim;

        using Base::get_nq;
        using Base::get_nq_dim;

        using Base::get_nvb;
        using Base::get_nvb_dim;

        using Base::get_nvj;
        using Base::get_nvj_dim;

        using Base::get_nv;
        using Base::get_nv_dim;

        using Base::get_nrotors;
        using Base::get_nrotors_dim;

        using Base::get_nx;
        using Base::get_nx_dim;

        using Base::get_ndx;
        using Base::get_ndx_dim;

        using Base::get_nua;
        using Base::get_nua_dim;

    protected:
        VectorNx_t lb_;
        VectorNx_t ub_;

    }; // class StateEuclideanTpl

} // namespace galileo

#endif // __galileo_core_states_euclidean_hpp__
