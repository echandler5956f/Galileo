#ifndef __galileo_core_node_collections_hpp__
#define __galileo_core_node_collections_hpp__

#include "galileo/core/node/fwd.hpp"
// #include "galileo/core/node/nodes.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename _Scalar, int _Options>
    struct NodeCollectionDefaultTpl
    {
        using Scalar = _Scalar;
        enum
        {
            Options = _Options
        };

        using NodeModelVariant = boost::variant<NodeModelVoid>;
        using NodeDataVariant = boost::variant<NodeDataVoid>;
    };

    using NodeModelVariant = typename NodeCollectionDefault::NodeModelVariant;
    using NodeDataVariant = typename NodeCollectionDefault::NodeDataVariant;

} // namespace galileo

#endif // __galileo_core_node_collections_hpp__
