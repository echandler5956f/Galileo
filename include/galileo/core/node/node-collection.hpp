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
        typedef _Scalar Scalar;
        enum
        {
            Options = _Options
        };

        typedef boost::variant<
            NodeModelVoid>
            NodeModelVariant;

        typedef boost::variant<
            NodeDataVoid>
            NodeDataVariant;
    };

    typedef NodeCollectionDefault::NodeModelVariant NodeModelVariant;
    typedef NodeCollectionDefault::NodeDataVariant NodeDataVariant;

} // namespace galileo

#endif // __galileo_core_node_collections_hpp__
