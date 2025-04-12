#ifndef __galileo_predictive_references_reference_manager_base_hpp__
#define __galileo_predictive_references_reference_manager_base_hpp__

#include "galileo/predictive/fwd.hpp"

#include "galileo/utils/threads/buffered_value.hpp"
#include "galileo/predictive/references/mode-schedule-base.hpp"
#include "galileo/predictive/references/target-trajectories-base.hpp"

#define GALILEO_REFERENCE_MANAGER_BASIC_TYPEDEF(ReferenceManager) \
    using VarScalar = typename traits<PrimalSolution>::VarScalar; \
    using NumScalar = typename traits<PrimalSolution>::NumScalar; \
    static constexpr int Options = traits<PrimalSolution>::Options;

#define GALILEO_REFERENCE_MANAGER_TYPEDEF(ReferenceManager)                               \
    using ModeSchedule_t = typename traits<ReferenceManager>::ModeSchedule_t;             \
    using TargetTrajectories_t = typename traits<ReferenceManager>::TargetTrajectories_t; \
    using VectorXv = typename traits<ReferenceManager>::VectorXv;                         \
    using VectorXn = typename traits<ReferenceManager>::VectorXn;

namespace galileo
{

    template <typename Derived>
    class ReferenceManagerBase : internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using ReferenceManagerDerived = typename traits<Derived>::ReferenceManagerDerived;
        GALILEO_REFERENCE_MANAGER_BASIC_TYPEDEF(ReferenceManagerDerived);
        GALILEO_REFERENCE_MANAGER_TYPEDEF(ReferenceManagerDerived);

        template <typename StateVectorType>
        void preSolverRun(NumScalar init_time, const Eigen::MatrixBase<StateVectorType> &init_state, NumScalar final_time)
        {
            target_trajectories_.updateFromBuffer();
            mode_schedule_.updateFromBuffer();
            modifyReferences(init_time, init_state, final_time, target_trajectories_.get(), mode_schedule_.get());
        }

        const ModeSchedule_t &getModeSchedule() const { return mode_schedule_.get(); }

        void setModeSchedule(const ModeSchedule_t &mode_schedule) { mode_schedule_.setBuffer(mode_schedule); }
        void setModeSchedule(ModeSchedule_t &&mode_schedule) { mode_schedule_.setBuffer(std::move(mode_schedule)); }

        const TargetTrajectories_t &getTargetTrajectories() const { return target_trajectories_.get(); }

        void setTargetTrajectories(const TargetTrajectories_t &target_trajectories) { target_trajectories_.setBuffer(target_trajectories); }
        void setTargetTrajectories(TargetTrajectories_t &&target_trajectories) { target_trajectories_.setBuffer(std::move(target_trajectories)); }

    protected:
        template <typename StateVectorType>
        void modifyReferences(NumScalar init_time, const Eigen::MatrixBase<StateVectorType> &init_state, NumScalar final_time, const TargetTrajectories_t &target_trajectories, const ModeSchedule_t &mode_schedule)
        {
        }

        inline ReferenceManagerBase()
        {
        }

        inline ReferenceManagerBase(const ReferenceManagerBase &clone)
        {
            *this = clone;
        }

        inline ReferenceManagerBase &operator=(const ReferenceManagerBase &clone)
        {
            return *this;
        }

        BufferedValue<ModeSchedule_t> mode_schedule_;
        BufferedValue<TargetTrajectories_t> target_trajectories_;

    }; // class ReferenceManagerBase

} // namespace galileo

#endif // __galileo_predictive_references_reference_manager_base_hpp__