#ifndef __galileo_multibody_predictive_nodes_fwd_hpp__
#define __galileo_multibody_predictive_nodes_fwd_hpp__

#include "galileo/domains/multibody/predictive/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct NodeModelContactFwdDynTpl;
    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct NodeDataContactFwdDynTpl;

} // namespace galileo

#endif // __galileo_multibody_predictive_nodes_fwd_hpp__
