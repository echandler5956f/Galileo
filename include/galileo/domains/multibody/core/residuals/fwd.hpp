#ifndef __galileo_multibody_core_residuals_fwd_hpp__
#define __galileo_multibody_core_residuals_fwd_hpp__

#include "galileo/domains/multibody/core/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    class ResidualModelFramePlacementTpl;
    template <typename PhaseSpec>
    struct ResidualDataFramePlacementTpl;

    template <typename PhaseSpec>
    class ResidualModelFrameTranslationTpl;
    template <typename PhaseSpec>
    struct ResidualDataFrameTranslationTpl;

    template <typename PhaseSpec>
    class ResidualModelFrameVelocityTpl;
    template <typename PhaseSpec>
    struct ResidualDataFrameVelocityTpl;

    template <typename PhaseSpec>
    class ResidualModelCoMPositionTpl;
    template <typename PhaseSpec>
    struct ResidualDataCoMPositionTpl;

    template <typename PhaseSpec>
    class ResidualModelMultibodyStateTpl;

} // namespace galileo

#endif // __galileo_multibody_core_residuals_fwd_hpp__
