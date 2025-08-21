#ifndef __galileo_core_states_state_euclidean_hpp__
#define __galileo_core_states_state_euclidean_hpp__

#include "galileo/core/states/state-base.hpp"

namespace galileo
{

    template <typename Spec>
    class StateEuclideanTpl : public StateBase<StateEuclideanTpl<Spec>, Spec>
    {
    public:
        using SS = Spec;
        GALILEO_SYSTEM_SPEC_MASTER_TYPEDEF(SS);

        using Base = StateBase<StateEuclideanTpl<SS>, SS>;

        StateEuclideanTpl(const SS &spec, const VectorNx_t &lb, const VectorNx_t &ub) : Base(spec), lb_(lb), ub_(ub)
        {
            GALILEO_ASSERT(IsValidEuclidean(*this),
                           "StateEuclideanTpl: Invalid dimensions for Euclidean state (nx != ndx)");
        }

        VectorNx_t zero() const { return VectorNx_t::Zero(get_nx()); }
        VectorNx_t rand() const { return VectorNx_t::Random(get_nx()); }

        template <typename StateVector1, typename StateVector2, typename StateTangentVector>
        void diff(const Eigen::MatrixBase<StateVector1> &x0,
                  const Eigen::MatrixBase<StateVector2> &x1,
                  Eigen::MatrixBase<StateTangentVector> &dxout) const
        {
            dxout = x1 - x0;
        }

        template <typename StateVector1, typename StateTangentVector, typename StateVector2>
        void integrate(const Eigen::MatrixBase<StateVector1> &x,
                       const Eigen::MatrixBase<StateTangentVector> &dx,
                       Eigen::MatrixBase<StateVector2> &xout) const
        {
            xout = x + dx;
        }

        template <Jcomponent jc = BOTH,
                  typename StateVector1,
                  typename StateVector2,
                  typename JMatrix1,
                  typename JMatrix2>
        void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                   const Eigen::MatrixBase<StateVector2> &x1,
                   Eigen::MatrixBase<JMatrix1> &Jfirst,
                   Eigen::MatrixBase<JMatrix2> &Jsecond) const
        {
            if constexpr (IsFirst<jc> || IsBoth<jc>)
            {
                Jfirst.setZero();
                Jfirst.diagonal() = VectorNdx_t::Constant(get_ndx(), VarScalar(-1.));
            }
            if constexpr (IsSecond<jc> || IsBoth<jc>)
            {
                Jsecond.setZero();
                Jsecond.diagonal() = VectorNdx_t::Constant(get_ndx(), VarScalar(1.));
            }
        }

        template <Jcomponent jc = BOTH,
                  AssignmentOp op = SETTO,
                  typename StateVector,
                  typename StateTangentVector,
                  typename JMatrix1,
                  typename JMatrix2>
        void Jintegrate(const Eigen::MatrixBase<StateVector> &x,
                        const Eigen::MatrixBase<StateTangentVector> &dx,
                        Eigen::MatrixBase<JMatrix1> &Jfirst,
                        Eigen::MatrixBase<JMatrix2> &Jsecond) const
        {
            if constexpr (IsFirst<jc> || IsBoth<jc>)
            {
                if constexpr (IsSetTo<op>)
                    Jfirst.diagonal().array() = VarScalar(1.);
                else if constexpr (IsAddTo<op>)
                    Jfirst.diagonal().array() += VarScalar(1.);
                else if constexpr (IsRmFrom<op>)
                    Jfirst.diagonal().array() -= VarScalar(1.);
            }

            if constexpr (IsSecond<jc> || IsBoth<jc>)
            {
                if constexpr (IsSetTo<op>)
                    Jsecond.diagonal().array() = VarScalar(1.);
                else if constexpr (IsAddTo<op>)
                    Jsecond.diagonal().array() += VarScalar(1.);
                else if constexpr (IsRmFrom<op>)
                    Jsecond.diagonal().array() -= VarScalar(1.);
            }
        }

        template <Jcomponent jc, typename StateVector, typename StateTangentVector, typename JMatrix>
        void JintegrateTransport(const Eigen::MatrixBase<StateVector> &x,
                                 const Eigen::MatrixBase<StateTangentVector> &dx,
                                 Eigen::MatrixBase<JMatrix> &Jin) const
        {
            // Nothing to do
        }

        const VectorNx_t &get_lb() const { return lb_; }
        const VectorNx_t &get_ub() const { return ub_; }

        template <typename StateVector>
        void set_lb(const Eigen::MatrixBase<StateVector> &lb)
        {
            lb_ = lb;
        }

        template <typename StateVector>
        void set_ub(const Eigen::MatrixBase<StateVector> &ub)
        {
            ub_ = ub;
        }

        using Base::diff_dx;
        using Base::integrate_x;
        using Base::Jdiff_Js;
        using Base::Jintegrate_Js;
        using Base::get_spec;
        using Base::get_nq;
        using Base::get_nq_dim;
        using Base::get_nv;
        using Base::get_nv_dim;
        using Base::get_nx;
        using Base::get_nx_dim;
        using Base::get_ndx;
        using Base::get_ndx_dim;
        using Base::get_nua;
        using Base::get_nua_dim;

        /**
         * @brief Display Euclidean state information to output stream
         */
        void display(std::ostream &os, const std::string &indent = "  ") const
        {
            os << indent << "Spec: {\n";
            get_spec().display(os, indent + "  ");
            os << indent << "}\n";
            os << indent << "State Bounds: {\n";
            os << indent << "  Lower bounds: " << lb_.transpose() << "\n";
            os << indent << "  Upper bounds: " << ub_.transpose() << "\n";
            os << indent << "}\n";
            os << indent << "Euclidean configuration: " << (IsValidEuclidean(*this) ? "VALID" : "INVALID") << "\n";
        }

    protected:
        VectorNx_t lb_;
        VectorNx_t ub_;

    }; // class StateEuclideanTpl

    template <typename SS>
    bool IsValidEuclidean(const StateEuclideanTpl<SS> &state)
    {
        return state.get_nx() == state.get_ndx();
    }

} // namespace galileo

#endif // __galileo_core_states_state_euclidean_hpp__
