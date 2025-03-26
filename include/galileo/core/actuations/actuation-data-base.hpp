#ifndef __galileo_core_actuations_actuation_data_base_hpp__
#define __galileo_core_actuations_actuation_data_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/core/actuations/actuation-model-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        struct ActuationDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActuationDerived = typename traits<Derived>::ActuationDerived;
            GALILEO_ACTUATIONS_BASIC_TYPEDEF(ActuationDerived);
            GALILEO_ACTUATIONS_CONSTANTS(ActuationDerived);
            GALILEO_ACTUATIONS_DATA_TYPEDEF(ActuationDerived);

        protected:
            inline ActuationDataBase()
            {
            }

            inline ActuationDataBase(const ActuationDataBase &clone)
            {
                *this = clone;
            }

            inline ActuationDataBase &operator=(const ActuationDataBase &clone)
            {
                return *this;
            }

        }; // struct ActuationDataBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_actuation_data_base_hpp__
