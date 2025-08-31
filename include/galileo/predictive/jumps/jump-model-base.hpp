#ifndef __galileo_predictive_jumps_jump_model_base_hpp__
#define __galileo_predictive_jumps_jump_model_base_hpp__

#include "galileo/predictive/jumps/jump-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class JumpModelBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using State_t = typename PS::State_t;
        using CostModelManager_t = typename PS::CostModelManager_t;
        using ConstraintModelManager_t = typename PS::ConstraintModelManager_t;

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x);
        }

        Data_t createData(MemoryArena &arena) const { return this->derived().createData(arena); }

        const PS &get_ps() const { return ps_; }
        const State_t &get_state() const { return state_; }
        const CostModelManager_t &get_costs() const { return this->derived().get_costs(); }
        const ConstraintModelManager_t &get_constraints() const { return this->derived().get_constraints(); }

    protected:
        inline JumpModelBase(const PS &ps, const State_t &state) : ps_(ps), state_(state) {}
        inline JumpModelBase(const JumpModelBase &clone) : ps_(clone.ps_), state_(clone.state_) {}
        inline JumpModelBase &operator=(const JumpModelBase &clone)
        {
            ps_ = clone.ps_;
            state_ = clone.state_;
            return *this;
        }

        PS ps_;
        State_t state_;

    }; // class JumpModelBase

} // namespace galileo

#endif // __galileo_predictive_jumps_jump_model_base_hpp__
