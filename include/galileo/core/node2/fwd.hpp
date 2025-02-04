#ifndef __galileo_core_node_fwd_hpp__
#define __galileo_core_node_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    struct NodeModelVoid
    {
    };

    struct NodeDataVoid
    {
    };

    template <typename Scalar, int Options = context::Options>
    struct NodeCollectionDefaultTpl;
    using NodeCollectionDefault = NodeCollectionDefaultTpl<context::Scalar>;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class NodeCollectionTpl = NodeCollectionDefaultTpl>
    struct NodeModelTpl;
    using NodeModel = NodeModelTpl<context::Scalar>;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class NodeCollectionTpl = NodeCollectionDefaultTpl>
    struct NodeDataTpl;
    using NodeData = NodeDataTpl<context::Scalar>;

} // namespace galileo

#include "galileo/core/fwd.hpp"

#endif // __galileo_core_node_fwd_hpp__
