#ifndef __galileo_core_residuals_residual_data_base_hpp__
#define __galileo_core_residuals_residual_data_base_hpp__

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/core/residuals/residual-model-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        struct ResidualDataBase : internal::CRTP<ResidualDataBase<Derived>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ResidualDerived = typename traits<Derived>::ResidualDerived;
            GALILEO_RESIDUAL_BASIC_TYPEDEF(ResidualDerived);
            GALILEO_RESIDUAL_CONSTANTS(ResidualDerived);
            GALILEO_RESIDUAL_DATA_TYPEDEF(ResidualDerived);

            // R_t R;
            // Rx_t Rx;
            // Ru_t Ru;
            // Arr_Rx_t Arr_Rx;
            // Arr_Ru_t Arr_Ru;

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

    } // namespace core

} // namespace galileo

#endif // __galileo_core_residuals_residual_data_base_hpp__
