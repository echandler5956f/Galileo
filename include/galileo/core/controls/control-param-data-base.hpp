#ifndef __galileo_core_controls_control_param_data_base_hpp__
#define __galileo_core_controls_control_param_data_base_hpp__

#include "galileo/core/controls/control-param-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ControlParamDataTpl
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = typename PS::ControlParamMeta_t;
        using Model_t = typename PS::ControlParamModel_t;
        using Data_t = typename PS::ControlParamData_t;

        ControlParamDataTpl(const Model_t &model)
            : u(model.get_nu()),
              w(model.get_nw()),
              du_dw(model.get_nu(),
                    model.get_nw())
        {
            u.setZero();
            w.setZero();
            du_dw.setZero();
        }

        U_t u;
        W_t w;
        Uw_t du_dw;

    }; // struct ControlParamDataTpl

} // namespace galileo

#endif // __galileo_core_controls_control_param_data_base_hpp__
