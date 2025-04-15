#ifndef __galileo_core_residuals_residual_data_base_hpp__
#define __galileo_core_residuals_residual_data_base_hpp__

#include "galileo/core/residuals/residual-base.hpp"

#define GALILEO_RESIDUAL_DATA_TYPEDEF(Residual)           \
    using R_t = typename traits<Residual>::R_t;           \
    using Rx_t = typename traits<Residual>::Rx_t;         \
    using Ru_t = typename traits<Residual>::Ru_t;         \
    using Arr_Rx_t = typename traits<Residual>::Arr_Rx_t; \
    using Arr_Ru_t = typename traits<Residual>::Arr_Ru_t;

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct ResidualDataBase : internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        FORWARD_ACCESSOR(R_t, R);
        FORWARD_ACCESSOR(Rx_t, Rx);
        FORWARD_ACCESSOR(Ru_t, Ru);
        FORWARD_ACCESSOR(Arr_Rx_t, Arr_Rx);
        FORWARD_ACCESSOR(Arr_Ru_t, Arr_Ru);

    protected:
        inline ResidualDataBase()
        {
        }

        inline ResidualDataBase(const ResidualDataBase &clone)
        {
            *this = clone;
        }

        inline ResidualDataBase &operator=(const ResidualDataBase &clone)
        {
            return *this;
        }

    }; // struct ResidualDataBase

} // namespace galileo

#endif // __galileo_core_residuals_residual_data_base_hpp__
