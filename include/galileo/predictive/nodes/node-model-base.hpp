#ifndef __galileo_predictive_nodes_node_model_base_hpp__
#define __galileo_predictive_nodes_node_model_base_hpp__

#include "galileo/predictive/nodes/node-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class NodeModelBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using DimNU_t = typename traits<Meta_t>::DimNU_t;

        using State_t = typename PS::State_t;
        using CostModelManager_t = typename PS::CostModelManager_t;
        using ConstraintModelManager_t = typename PS::ConstraintModelManager_t;

        using NumScalar = typename PS::NumScalar;
        using UBound_t = Eigen::GMatrix<NumScalar, DimNU_t::Value, 1, PS::Options>;

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calc(data, x, u);
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x, u);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void quasiStatic(Data_t &data,
                         const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlVectorType> &u,
                         const int maxiter,
                         const NumScalar tol) const
        {
            this->derived().quasiStatic(data, x, u, maxiter, tol);
        }

        Data_t createData(MemoryArena &arena) const { return this->derived().createData(arena); }

        const PS &get_ps() const { return ps_; }
        const State_t &get_state() const { return state_; }
        const CostModelManager_t &get_costs() const { return this->derived().get_costs(); }
        const ConstraintModelManager_t &get_constraints() const { return this->derived().get_constraints(); }
        const UBound_t &get_u_lb() const { return u_lb_; }
        const UBound_t &get_u_ub() const { return u_ub_; }
        void set_u_lb(const UBound_t &u_lb) { u_lb_ = u_lb; }
        void set_u_ub(const UBound_t &u_ub) { u_ub_ = u_ub; }

    protected:
        inline NodeModelBase(PS &ps, const State_t &state)
            : ps_(ps),
              state_(state),
              u_lb_(UBound_t::Zero(get_ps().get_nu())),
              u_ub_(UBound_t::Zero(get_ps().get_nu()))
        {
        }
        inline NodeModelBase(const NodeModelBase &clone)
            : ps_(clone.ps_), state_(clone.state_), u_lb_(clone.u_lb_), u_ub_(clone.u_ub_)
        {
        }
        inline NodeModelBase &operator=(const NodeModelBase &clone)
        {
            ps_ = clone.ps_;
            state_ = clone.state_;
            u_lb_ = clone.u_lb_;
            u_ub_ = clone.u_ub_;
            return *this;
        }

        PS ps_;
        State_t state_;
        UBound_t u_lb_;
        UBound_t u_ub_;

    }; // class NodeModelBase

} // namespace galileo

#endif // __galileo_predictive_nodes_node_model_base_hpp__
