#ifndef __galileo_core_node_basic_visitors_hxx__
#define __galileo_core_node_basic_visitors_hxx__

#include "galileo/core/node/node-basic-visitors.hpp"
#include "galileo/core/visitor.hpp"

namespace galileo
{

  template<typename Scalar, int Options, template<typename S, int O> class NodeCollectionTpl>
  struct CreateNodeData : boost::static_visitor<NodeDataTpl<Scalar, Options, NodeCollectionTpl>>
  {
    typedef NodeCollectionTpl<Scalar, Options> NodeCollection;
    typedef typename NodeCollection::NodeModelVariant NodeModelVariant;
    typedef NodeDataTpl<Scalar, Options, NodeCollectionTpl> NodeDataVariant;

    template<typename NodeModelDerived>
    NodeDataVariant operator()(const NodeModelBase<NodeModelDerived> & node_model) const
    {
      return NodeDataVariant(node_model.createData());
    }

    static NodeDataVariant run(const NodeModelVariant & node_model)
    {
      return boost::apply_visitor(CreateNodeData(), node_model);
    }
  };

  template<typename Scalar, int Options, template<typename S, int O> class NodeCollectionTpl>
  inline NodeDataTpl<Scalar, Options, NodeCollectionTpl>
  createData(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> & node_model)
  {
    return CreateNodeData<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

} // namespace galileo

#endif // __galileo_core_node_basic_visitors_hxx__
