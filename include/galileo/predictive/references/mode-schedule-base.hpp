#ifndef __galileo_predictive_references_mode_schedule_base_hpp__
#define __galileo_predictive_references_mode_schedule_base_hpp__

#include "galileo/predictive/fwd.hpp"

#define GALILEO_MODE_SCHEDULE_BASIC_TYPEDEF(ModeSchedule) \
    using VarScalar = typename traits<ModeSchedule>::VarScalar;    \
    using NumScalar = typename traits<ModeSchedule>::NumScalar;    \
    static constexpr int Options = traits<ModeSchedule>::Options;

#define GALILEO_MODE_SCHEDULE_TYPEDEF(ModeSchedule)    \
    using VectorXv = typename traits<ModeSchedule>::VectorXv; \
    using VectorXn = typename traits<ModeSchedule>::VectorXn;

namespace galileo
{

    namespace predictive
    {

        template <typename Derived>
        class ModeScheduleBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ModeScheduleDerived = typename traits<Derived>::ModeScheduleDerived;
            GALILEO_MODE_SCHEDULE_BASIC_TYPEDEF(ModeScheduleDerived);
            GALILEO_MODE_SCHEDULE_TYPEDEF(ModeScheduleDerived);

            size_t getMode(NumScalar t) const
            {
                return derived().getMode(t);
            }

            std::vector<NumScalar> event_times_;
            std::vector<size_t> mode_sequence_;

        }; // class ModeScheduleBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_references_mode_schedule_base_hpp__
