#ifndef __galileo_predictive_jumps_jump_model_base_hpp__
#define __galileo_predictive_jumps_jump_model_base_hpp__

#include "galileo/predictive/jumps/jump-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class JumpModelBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

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

        Data_t createData() const { return this->derived().createData(); }

        const PS &get_ps() const { return ps_.get(); }
        const State_t &get_state() const { return get_ps().get_state(); }
        const RobotModel_t &get_robot() const { return robot_.get(); }
        const CostModelManager_t &get_costs() const { return this->derived().get_costs(); }
        const ConstraintModelManager_t &get_constraints() const { return this->derived().get_constraints(); }

    protected:
        inline JumpModelBase(const PS &ps) : ps_(ps), robot_(ps.get_state().get_robot()) {}
        inline JumpModelBase(const JumpModelBase &clone) : ps_(clone.ps_), robot_(clone.robot_) {}
        inline JumpModelBase &operator=(const JumpModelBase &clone)
        {
            ps_ = clone.ps_;
            robot_ = clone.robot_;
            return *this;
        }

        std::reference_wrapper<const PS> ps_;
        std::reference_wrapper<const RobotModel_t> robot_;

    }; // class JumpModelBase

} // namespace galileo

#endif // __galileo_predictive_jumps_jump_model_base_hpp__
