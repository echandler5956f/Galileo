#ifndef __galileo_core_controls_fwd_hpp__
#define __galileo_core_controls_fwd_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ControlParamDataTpl;

    template <typename PhaseSpec, int NOrder_>
    struct ControlParamModelJacobiPolynomialTpl;

} // namespace galileo

#endif // __galileo_core_controls_fwd_hpp__