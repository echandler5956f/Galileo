#ifndef __galileo_predictive_phases_phase_data_base_hpp__
#define __galileo_predictive_phases_phase_data_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct PhaseDataBase
        : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

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

} // namespace galileo

#endif // __galileo_predictive_phases_phase_data_base_hpp__
