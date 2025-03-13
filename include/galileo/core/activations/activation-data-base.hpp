#ifndef __galileo_core_activations_activation_data_base_hpp__
#define __galileo_core_activations_activation_data_base_hpp__

#include "galileo/core/activations/activation-base.hpp"
#include "galileo/core/activations/activation-model-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        struct ActivationDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActivationDerived = typename traits<Derived>::ActivationDerived;
            GALILEO_ACTIVATION_BASIC_TYPEDEF(ActivationDerived);
            GALILEO_ACTIVATION_CONSTANTS(ActivationDerived);
            GALILEO_ACTIVATION_DATA_TYPEDEF(ActivationDerived);

            A_t A;
            Ar_t Ar;
            Arr_t Arr;

        protected:
            inline ActivationDataBase()
            {
            }

            inline ActivationDataBase(const ActivationDataBase &clone)
            {
                *this = clone;
            }

            inline ActivationDataBase &operator=(const ActivationDataBase &clone)
            {
                return *this;
            }

        }; // struct ActivationDataBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_activations_activation_data_base_hpp__
