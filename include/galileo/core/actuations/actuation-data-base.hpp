#ifndef __galileo_core_actuations_actuation_data_base_hpp__
#define __galileo_core_actuations_actuation_data_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

#define GALILEO_ACTUATION_DATA_TYPEDEF(Actuation) \
    using VectorNv_t = typename traits<Actuation>::VectorNv_t; \
    using VectorNua_t = typename traits<Actuation>::VectorNua_t; \
    using MatrixNvNdx_t = typename traits<Actuation>::MatrixNvNdx_t; \
    using MatrixNvNua_t = typename traits<Actuation>::MatrixNvNua_t; \
    using MatrixNuaNv_t = typename traits<Actuation>::MatrixNuaNv_t; \
    using BoolArrayNv_t = typename traits<Actuation>::BoolArrayNv_t;

namespace galileo
{

    template <typename Derived, typename SystemSpec>
    struct ActuationDataBase : public internal::CRTP<Derived>
    {
    public:
        using SS = SystemSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_ACTUATION_DATA_TYPEDEF(Meta_t);

        FORWARD_ACCESSOR(VectorNv_t, tau);
        FORWARD_ACCESSOR(VectorNua_t, u);
        FORWARD_ACCESSOR(MatrixNvNdx_t, dtau_dx);
        FORWARD_ACCESSOR(MatrixNvNua_t, dtau_du);
        FORWARD_ACCESSOR(MatrixNuaNv_t, Mtau);
        FORWARD_ACCESSOR(BoolArrayNv_t, tau_set);

    protected:
        inline ActuationDataBase() {}
        inline ActuationDataBase(const ActuationDataBase &clone) { *this = clone; }
        inline ActuationDataBase &operator=(const ActuationDataBase &clone) { return *this; }

    }; // struct ActuationDataTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_data_base_hpp__
