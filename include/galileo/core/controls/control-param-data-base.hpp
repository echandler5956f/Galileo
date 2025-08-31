#ifndef __galileo_core_controls_control_param_data_base_hpp__
#define __galileo_core_controls_control_param_data_base_hpp__

#include "galileo/core/controls/control-param-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct ControlParamDataBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using VectorNu_t = ArenaMatrixTpl<typename PS::VectorNu_t>;
        using VectorNw_t = ArenaMatrixTpl<typename PS::VectorNw_t>;
        using MatrixNuNw_t = ArenaMatrixTpl<typename PS::MatrixNuNw_t>;

        FORWARD_ACCESSOR(VectorNu_t, u);
        FORWARD_ACCESSOR(VectorNw_t, w);
        FORWARD_ACCESSOR(MatrixNuNw_t, du_dw);

    protected:
        inline ControlParamDataBase() {}
        inline ControlParamDataBase(const ControlParamDataBase &clone) { *this = clone; }
        inline ControlParamDataBase &operator=(const ControlParamDataBase &clone) { return *this; }

    }; // struct ControlParamDataTpl

} // namespace galileo

#endif // __galileo_core_controls_control_param_data_base_hpp__
