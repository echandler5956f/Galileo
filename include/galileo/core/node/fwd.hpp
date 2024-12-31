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
    typedef NodeCollectionDefaultTpl<context::Scalar> NodeCollectionDefault;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class NodeCollectionTpl = NodeCollectionDefaultTpl>
    struct NodeModelTpl;
    typedef NodeModelTpl<context::Scalar> NodeModel;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class NodeCollectionTpl = NodeCollectionDefaultTpl>
    struct NodeDataTpl;
    typedef NodeDataTpl<context::Scalar> NodeData;

} // namespace galileo

#include "galileo/core/fwd.hpp"

#endif // __galileo_core_node_fwd_hpp__
