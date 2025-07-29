#ifndef __galileo_predictive_nodes_fwd_hpp__
#define __galileo_predictive_nodes_fwd_hpp__

#include "galileo/predictive/fwd.hpp"

namespace galileo
{

    template <
        typename PhaseSpec,
        template <typename PS> class ContactCollectionTpl>
    struct NodeModelContactFwdDynTpl;
    template <
        typename PhaseSpec,
        template <typename PS> class ContactCollectionTpl>
    struct NodeDataContactFwdDynTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class ImpulseCollectionTpl>
    struct NodeModelImpulseFwdDynTpl;
    template <
        typename PhaseSpec,
        template <typename PS> class ImpulseCollectionTpl>
    struct NodeDataImpulseFwdDynTpl;

} // namespace galileo

#endif // __galileo_predictive_nodes_fwd_hpp__
