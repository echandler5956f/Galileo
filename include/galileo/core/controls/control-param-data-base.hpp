#ifndef __galileo_core_controls_control_param_data_base_hpp__
#define __galileo_core_controls_control_param_data_base_hpp__

#include "galileo/core/controls/control-param-base.hpp"
#include "galileo/core/controls/control-param-model-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ControlParamDataTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename PS::ControlParamMeta_t;
        using Model_t = typename PS::ControlParamModel_t;
        using Data_t = typename PS::ControlParamData_t;

        typename PS::U_t u;
        typename PS::W_t w;
        typename PS::Uw_t du_dw;

    }; // struct ControlParamDataTpl

} // namespace galileo

#endif // __galileo_core_controls_control_param_data_base_hpp__
