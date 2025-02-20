#ifndef __galileo_predictive_phases_phase_data_base_hpp__
#define __galileo_predictive_phases_phase_data_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-model-base.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename Derived>
        struct PhaseDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using NodeDerived = typename traits<Derived>::NodeDerived;
            using SegmentDerived = typename traits<Derived>::SegmentDerived;
            using PhaseDerived = typename traits<Derived>::PhaseDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_BASIC_TYPEDEF(SegmentDerived);
            GALILEO_PHASE_BASIC_TYPEDEF(PhaseDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_CONSTANTS(SegmentDerived);
            GALILEO_PHASE_CONSTANTS(PhaseDerived);

            GALILEO_NODE_DATA_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_DATA_TYPEDEF(SegmentDerived);
            GALILEO_PHASE_DATA_TYPEDEF(PhaseDerived);

            SegmentDataVector segments;
            JumpData_t jump;

        protected:
            inline PhaseDataBase()
            {
            }

            inline PhaseDataBase(const PhaseDataBase &clone)
            {
                *this = clone;
            }

            inline PhaseDataBase &operator=(const PhaseDataBase &clone)
            {
                return *this;
            }

        }; // struct PhaseDataBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_data_base_hpp__
