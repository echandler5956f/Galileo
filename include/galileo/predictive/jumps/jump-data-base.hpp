#ifndef __galileo_predictive_jumps_jump_data_base_hpp__
#define __galileo_predictive_jumps_jump_data_base_hpp__

#include "galileo/predictive/jumps/jump-base.hpp"
#include "galileo/predictive/jumps/jump-model-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct JumpDataBase
        : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        FORWARD_ACCESSOR(CostDataManager_t, costs);
        FORWARD_ACCESSOR(ConstraintDataManager_t, constraints);

        FORWARD_ACCESSOR(XNext_t, XNext);
        FORWARD_ACCESSOR(XNextx_t, XNextx);

        FORWARD_ACCESSOR(L_t, L);
        FORWARD_ACCESSOR(Lx_t, Lx);
        FORWARD_ACCESSOR(Lxx_t, Lxx);

        FORWARD_ACCESSOR(H_t, H);
        FORWARD_ACCESSOR(Hx_t, Hx);

        FORWARD_ACCESSOR(G_t, G);
        FORWARD_ACCESSOR(Gx_t, Gx);

    protected:
        inline JumpDataBase()
        {
        }

        inline JumpDataBase(const JumpDataBase &clone)
        {
            *this = clone;
        }

        inline JumpDataBase &operator=(const JumpDataBase &clone)
        {
            return *this;
        }

    }; // struct JumpDataBase

} // namespace galileo

#endif // __galileo_predictive_jumps_jump_data_base_hpp__
