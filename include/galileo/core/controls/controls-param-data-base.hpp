#ifndef __galileo_core_controls_controls_param_data_base_hpp__
#define __galileo_core_controls_controls_param_data_base_hpp__

#include "galileo/core/controls/controls-param-base.hpp"
#include "galileo/core/controls/controls-param-model-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename PhaseSpec>
        struct ControlParamDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

        protected:
            inline ControlParamDataBase()
            {
            }

            inline ControlParamDataBase(const ControlParamDataBase &clone)
            {
                *this = clone;
            }

            inline ControlParamDataBase &operator=(const ControlParamDataBase &clone)
            {
                return *this;
            }

        }; // struct ControlParamDataBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_controls_controls_param_data_base_hpp__
