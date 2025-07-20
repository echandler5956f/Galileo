#ifndef __galileo_multibody_residuals_fwd_hpp__
#define __galileo_multibody_residuals_fwd_hpp__

#include "galileo/multibody/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    class ResidualModelControlTpl;
    template <typename PhaseSpec>
    struct ResidualDataControlTpl;

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
    class ResidualModelStateTpl;
    template <typename PhaseSpec>
    struct ResidualDataStateTpl;

} // namespace galileo

#endif // __galileo_multibody_residuals_fwd_hpp__
