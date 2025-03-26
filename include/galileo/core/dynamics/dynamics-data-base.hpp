#ifndef __galileo_core_dynamics_dynamics_data_base_hpp__
#define __galileo_core_dynamics_dynamics_data_base_hpp__

#include "galileo/core/dynamics/dynamics-base.hpp"
#include "galileo/core/dynamics/dynamics-model-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        struct DynamicsDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using DynamicsDerived = typename traits<Derived>::DynamicsDerived;
            GALILEO_DYNAMICS_BASIC_TYPEDEF(DynamicsDerived);
            GALILEO_DYNAMICS_CONSTANTS(DynamicsDerived);
            GALILEO_DYNAMICS_DATA_TYPEDEF(DynamicsDerived);

            ActuationData_t actuation;
            F_t F;
            Fx_t Fx;
            Fu_t Fu;

        protected:
            inline DynamicsDataBase()
            {
            }

            inline DynamicsDataBase(const DynamicsDataBase &clone)
            {
                *this = clone;
            }

            inline DynamicsDataBase &operator=(const DynamicsDataBase &clone)
            {
                return *this;
            }

        }; // struct DynamicsDataBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_dynamics_dynamics_data_base_hpp__
