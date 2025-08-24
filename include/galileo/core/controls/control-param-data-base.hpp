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

        using Meta_t = typename PS::ControlParamMeta_t;
        using Model_t = typename PS::ControlParamModel_t;
        using Data_t = typename PS::ControlParamData_t;

        using U_t = typename PS::VectorNu_t;
        using W_t = typename PS::VectorNw_t;
        using Uw_t = typename PS::MatrixNuNw_t;

        ControlParamDataTpl(const Model_t &model)
            : u(model.get_ps().get_nu()),
              w(model.get_ps().get_nw()),
              du_dw(model.get_ps().get_nu(), model.get_ps().get_nw())
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
