#ifndef __galileo_core_actuations_actuation_data_base_hpp__
#define __galileo_core_actuations_actuation_data_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/core/actuations/actuation-model-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        struct ActuationsDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActuationsDerived = typename traits<Derived>::ActuationsDerived;
            GALILEO_ACTUATIONS_BASIC_TYPEDEF(ActuationsDerived);
            GALILEO_ACTUATIONS_CONSTANTS(ActuationsDerived);
            GALILEO_ACTUATIONS_DATA_TYPEDEF(ActuationsDerived);

        protected:
            inline ActuationsDataBase()
            {
            }

            inline ActuationsDataBase(const ActuationsDataBase &clone)
            {
                *this = clone;
            }

            inline ActuationsDataBase &operator=(const ActuationsDataBase &clone)
            {
                return *this;
            }

        }; // struct ActuationsDataBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_actuation_data_base_hpp__
