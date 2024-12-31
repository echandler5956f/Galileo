#ifndef __galileo_core_node_basic_visitors_hpp__
#define __galileo_core_node_basic_visitors_hpp__

#include "galileo/core/node/fwd.hpp"

namespace galileo
{

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline NodeDataTpl<Scalar, Options, NodeCollectionTpl>
    createData(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

} // namespace galileo

/* --- Details -------------------------------------------------------------------- */
// Included later
// #include "galileo/core/node/core-basic-visitors.hxx"

#endif // __galileo_core_node_basic_visitors_hpp__
