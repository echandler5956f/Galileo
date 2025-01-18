#ifndef __galileo_core_node_basic_visitors_hxx__
#define __galileo_core_node_basic_visitors_hxx__

#include "galileo/core/node/node-basic-visitors.hpp"
#include "galileo/core/visitor.hpp"

namespace galileo
{

  // Visitors on Node Models

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct CreateNodeData : boost::static_visitor<NodeDataTpl<Scalar, Options, NodeCollectionTpl>>
  {
    using NodeCollection = NodeCollectionTpl<Scalar, Options>;
    using NodeModelVariant = NodeCollection::NodeModelVariant;
    using NodeDataVariant = NodeDataTpl<Scalar, Options, NodeCollectionTpl>;

    template <typename NodeModelDerived>
    NodeDataVariant operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return NodeDataVariant(node_model.createData());
    }

    static NodeDataVariant run(const NodeModelVariant &node_model)
    {
      return boost::apply_visitor(CreateNodeData(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline NodeDataTpl<Scalar, Options, NodeCollectionTpl>
  createData(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return CreateNodeData<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename StateVectorType, typename ControlVectorType>
  struct NodeCalcVisitor
      : fusion::NodeUnaryVisitorBase<NodeCalcVisitor<StateVectorType, ControlVectorType>>
  {
    using ArgsType = boost::fusion::vector<const StateVectorType &, const ControlVectorType &>;

    template <typename NodeModel>
    static void algo(
        const galileo::NodeModelBase<NodeModel> &node_model,
        galileo::NodeDataBase<typename NodeModel::NodeDataDerived> &node_data,
        const Eigen::MatrixBase<StateVectorType> &xs,
        const Eigen::MatrixBase<ControlVectorType> &us)
    {
      node_model.calc(node_data.derived(), xs.derived(), us.derived());
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateVectorType, typename ControlVectorType>
  inline void calc(
      const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
      NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data,
      const Eigen::MatrixBase<StateVectorType> &xs,
      const Eigen::MatrixBase<ControlVectorType> &us)
  {
    using Algo = NodeCalcVisitor<StateVectorType, ControlVectorType>;

    Algo::run(node_model, node_data, typename Algo::ArgsType(xs.derived(), us.derived()));
  }

  template <typename StateVectorType, typename ControlVectorType>
  struct NodeCalcDiffVisitor
      : fusion::NodeUnaryVisitorBase<NodeCalcDiffVisitor<StateVectorType, ControlVectorType>>
  {
    using ArgsType = boost::fusion::vector<const StateVectorType &, const ControlVectorType &>;

    template <typename NodeModel>
    static void algo(
        const galileo::NodeModelBase<NodeModel> &node_model,
        galileo::NodeDataBase<typename NodeModel::NodeDataDerived> &node_data,
        const Eigen::MatrixBase<StateVectorType> &xs,
        const Eigen::MatrixBase<ControlVectorType> &us)
    {
      node_model.calcDiff(node_data.derived(), xs.derived(), us.derived());
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateVectorType, typename ControlVectorType>
  inline void calcDiff(
      const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
      NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data,
      const Eigen::MatrixBase<StateVectorType> &xs,
      const Eigen::MatrixBase<ControlVectorType> &us)
  {
    using Algo = NodeCalcDiffVisitor<StateVectorType, ControlVectorType>;

    Algo::run(node_model, node_data, typename Algo::ArgsType(xs.derived(), us.derived()));
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetXLBVisitor
      : boost::static_visitor<typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::X_t>
  {
    using ReturnType = typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::X_t;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_x_lb();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetXLBVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::X_t
  node_get_x_lb(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetXLBVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetXUBVisitor
      : boost::static_visitor<typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::X_t>
  {
    using ReturnType = typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::X_t;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_x_ub();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetXUBVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::X_t
  node_get_x_ub(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetXUBVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetULBVisitor
      : boost::static_visitor<typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::U_t>
  {
    using ReturnType = typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::U_t;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_u_lb();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetULBVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::U_t
  node_get_u_lb(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetULBVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetUUBVisitor
      : boost::static_visitor<typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::U_t>
  {
    using ReturnType = typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::U_t;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_u_ub();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetUUBVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::U_t
  node_get_u_ub(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetUUBVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetHLBVisitor
      : boost::static_visitor<typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::H_t>
  {
    using ReturnType = typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::H_t;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_h_lb();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetHLBVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::H_t
  node_get_h_lb(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetHLBVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetHUBVisitor
      : boost::static_visitor<typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::H_t>
  {
    using ReturnType = typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::H_t;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_h_ub();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetHUBVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::H_t
  node_get_h_ub(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetHUBVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetGLBVisitor
      : boost::static_visitor<typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::G_t>
  {
    using ReturnType = typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::G_t;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_g_lb();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetGLBVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::G_t
  node_get_g_lb(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetGLBVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetGUBVisitor
      : boost::static_visitor<typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::G_t>
  {
    using ReturnType = typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::G_t;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_g_ub();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetGUBVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::G_t
  node_get_g_ub(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetGUBVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetNHVisitor
      : boost::static_visitor<Eigen::Index>
  {
    using ReturnType = Eigen::Index;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_nh();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetNHVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline Eigen::Index
  node_get_nh(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetNHVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetNGVisitor
      : boost::static_visitor<Eigen::Index>
  {
    using ReturnType = Eigen::Index;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_ng();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetNGVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline Eigen::Index
  node_get_ng(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetNGVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGetIDVisitor
      : boost::static_visitor<Eigen::Index>
  {
    using ReturnType = Eigen::Index;

    template <typename NodeModelDerived>
    ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
    {
      return node_model.get_id();
    }

    static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
    {
      return boost::apply_visitor(NodeGetIDVisitor(), node_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline Eigen::Index
  node_get_id(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    return NodeGetIDVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateBoundVectorType>
  struct NodeSetXLBVisitor
      : boost::static_visitor<void>
  {
    using StateBoundVectorBase = Eigen::MatrixBase<StateBoundVectorType>;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, const StateBoundVectorBase &x_lb) const
    {
      node_model.set_x_lb(x_lb);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, const StateBoundVectorBase &x_lb)
    {
      boost::apply_visitor(NodeSetXLBVisitor(), node_model, x_lb);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateBoundVectorType>
  inline void
  node_set_x_lb(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                const Eigen::MatrixBase<StateBoundVectorType> &x_lb)
  {
    NodeSetXLBVisitor<Scalar, Options, NodeCollectionTpl, StateBoundVectorType>::run(node_model, x_lb.derived());
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateBoundVectorType>
  struct NodeSetXUBVisitor
      : boost::static_visitor<void>
  {
    using StateBoundVectorBase = Eigen::MatrixBase<StateBoundVectorType>;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, const StateBoundVectorBase &x_ub) const
    {
      node_model.set_x_ub(x_ub);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, const StateBoundVectorBase &x_ub)
    {
      boost::apply_visitor(NodeSetXUBVisitor(), node_model, x_ub);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateBoundVectorType>
  inline void
  node_set_x_ub(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                const Eigen::MatrixBase<StateBoundVectorType> &x_ub)
  {
    NodeSetXUBVisitor<Scalar, Options, NodeCollectionTpl, StateBoundVectorType>::run(node_model, x_ub.derived());
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename ControlBoundVectorType>
  struct NodeSetULBVisitor
      : boost::static_visitor<void>
  {
    using ControlBoundVectorBase = Eigen::MatrixBase<ControlBoundVectorType>;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, const ControlBoundVectorBase &u_lb) const
    {
      node_model.set_u_lb(u_lb);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, const ControlBoundVectorBase &u_lb)
    {
      boost::apply_visitor(NodeSetULBVisitor(), node_model, u_lb);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename ControlBoundVectorType>
  inline void
  node_set_u_lb(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                const Eigen::MatrixBase<ControlBoundVectorType> &u_lb)
  {
    NodeSetULBVisitor<Scalar, Options, NodeCollectionTpl, ControlBoundVectorType>::run(node_model, u_lb.derived());
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename ControlBoundVectorType>
  struct NodeSetUUBVisitor
      : boost::static_visitor<void>
  {
    using ControlBoundVectorBase = Eigen::MatrixBase<ControlBoundVectorType>;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, const ControlBoundVectorBase &u_ub) const
    {
      node_model.set_u_ub(u_ub);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, const ControlBoundVectorBase &u_ub)
    {
      boost::apply_visitor(NodeSetUUBVisitor(), node_model, u_ub);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename ControlBoundVectorType>
  inline void
  node_set_u_ub(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                const Eigen::MatrixBase<ControlBoundVectorType> &u_ub)
  {
    NodeSetUUBVisitor<Scalar, Options, NodeCollectionTpl, ControlBoundVectorType>::run(node_model, u_ub.derived());
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename EqualityBoundVectorType>
  struct NodeSetHLBVisitor
      : boost::static_visitor<void>
  {
    using EqualityBoundVectorBase = Eigen::MatrixBase<EqualityBoundVectorType>;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, const EqualityBoundVectorBase &h_lb) const
    {
      node_model.set_h_lb(h_lb);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, const EqualityBoundVectorBase &h_lb)
    {
      boost::apply_visitor(NodeSetHLBVisitor(), node_model, h_lb);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename EqualityBoundVectorType>
  inline void
  node_set_h_lb(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                const Eigen::MatrixBase<EqualityBoundVectorType> &h_lb)
  {
    NodeSetHLBVisitor<Scalar, Options, NodeCollectionTpl, EqualityBoundVectorType>::run(node_model, h_lb.derived());
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename EqualityBoundVectorType>
  struct NodeSetHUBVisitor
      : boost::static_visitor<void>
  {
    using EqualityBoundVectorBase = Eigen::MatrixBase<EqualityBoundVectorType>;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, const EqualityBoundVectorBase &h_ub) const
    {
      node_model.set_h_ub(h_ub);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, const EqualityBoundVectorBase &h_ub)
    {
      boost::apply_visitor(NodeSetHUBVisitor(), node_model, h_ub);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename EqualityBoundVectorType>
  inline void
  node_set_h_ub(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                const Eigen::MatrixBase<EqualityBoundVectorType> &h_ub)
  {
    NodeSetHUBVisitor<Scalar, Options, NodeCollectionTpl, EqualityBoundVectorType>::run(node_model, h_ub.derived());
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename InequalityBoundVectorType>
  struct NodeSetGLBVisitor
      : boost::static_visitor<void>
  {
    using InequalityBoundVectorBase = Eigen::MatrixBase<InequalityBoundVectorType>;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, const InequalityBoundVectorBase &g_lb) const
    {
      node_model.set_g_lb(g_lb);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, const InequalityBoundVectorBase &g_lb)
    {
      boost::apply_visitor(NodeSetGLBVisitor(), node_model, g_lb);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename InequalityBoundVectorType>
  inline void
  node_set_g_lb(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                const Eigen::MatrixBase<InequalityBoundVectorType> &g_lb)
  {
    NodeSetGLBVisitor<Scalar, Options, NodeCollectionTpl, InequalityBoundVectorType>::run(node_model, g_lb.derived());
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename InequalityBoundVectorType>
  struct NodeSetGUBVisitor
      : boost::static_visitor<void>
  {
    using InequalityBoundVectorBase = Eigen::MatrixBase<InequalityBoundVectorType>;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, const InequalityBoundVectorBase &g_ub) const
    {
      node_model.set_g_ub(g_ub);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, const InequalityBoundVectorBase &g_ub)
    {
      boost::apply_visitor(NodeSetGUBVisitor(), node_model, g_ub);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename InequalityBoundVectorType>
  inline void
  node_set_g_ub(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                const Eigen::MatrixBase<InequalityBoundVectorType> &g_ub)
  {
    NodeSetGUBVisitor<Scalar, Options, NodeCollectionTpl, InequalityBoundVectorType>::run(node_model, g_ub.derived());
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeSetNHVisitor
      : boost::static_visitor<void>
  {
    using ReturnType = Eigen::Index;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, Eigen::Index nh) const
    {
      node_model.set_nh(nh);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, Eigen::Index nh)
    {
      boost::apply_visitor(NodeSetNHVisitor(), node_model, nh);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline void
  node_set_nh(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
              Eigen::Index nh)
  {
    NodeSetNHVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model, nh);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeSetNGVisitor
      : boost::static_visitor<void>
  {
    using ReturnType = Eigen::Index;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, Eigen::Index ng) const
    {
      node_model.set_ng(ng);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, Eigen::Index ng)
    {
      boost::apply_visitor(NodeSetNGVisitor(), node_model, ng);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline void
  node_set_ng(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
              Eigen::Index ng)
  {
    NodeSetNGVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model, ng);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeSetIDVisitor
      : boost::static_visitor<void>
  {
    using ReturnType = Eigen::Index;

    template <typename NodeModelDerived>
    void operator()(NodeModelBase<NodeModelDerived> &node_model, Eigen::Index id) const
    {
      node_model.set_id(id);
    }

    static void run(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, Eigen::Index id)
    {
      boost::apply_visitor(NodeSetIDVisitor(), node_model, id);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline void
  node_set_id(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
              Eigen::Index id)
  {
    NodeSetIDVisitor<Scalar, Options, NodeCollectionTpl>::run(node_model, id);
  }

  template <typename NewScalar, typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeCastVisitor
      : fusion::NodeUnaryVisitorBase<
            NodeCastVisitor<NewScalar, Scalar, Options, NodeCollectionTpl>,
            typename CastType<NewScalar, NodeModelTpl<Scalar, Options, NodeCollectionTpl>>::type>
  {
    using ArgsType = fusion::NoArg;

    using ReturnType = typename CastType<NewScalar, NodeModelTpl<Scalar, Options, NodeCollectionTpl>>::type;

    template <typename NodeModelDerived>
    static ReturnType algo(const NodeModelBase<NodeModelDerived> &node_model)
    {
      return ReturnType(node_model.template cast<NewScalar>());
    }
  };

  template <typename NewScalar, typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  typename CastType<NewScalar, NodeModelTpl<Scalar, Options, NodeCollectionTpl>>::type
  cast_node(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
  {
    using Algo = NodeCastVisitor<NewScalar, Scalar, Options, NodeCollectionTpl>;
    return Algo::run(node_model);
  }

  // Visitors on Node Data

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeXDotVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::dX_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::dX_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_xdot();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeXDotVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::dX_t
  node_xdot(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeXDotVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeFxVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Fx_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Fx_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_fx();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeFxVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Fx_t
  node_fx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeFxVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeFuVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Fu_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Fu_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_fu();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeFuVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Fu_t
  node_fu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeFuVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeLVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::L_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::L_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_l();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeLVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::L_t
  node_l(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeLVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeLxVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lx_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lx_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_lx();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeLxVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lx_t
  node_lx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeLxVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeLuVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lu_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lu_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_lu();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeLuVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lu_t
  node_lu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeLuVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeLxxVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lxx_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lxx_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_lxx();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeLxxVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lxx_t
  node_lxx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeLxxVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeLxuVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lxu_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lxu_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_lxu();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeLxuVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lxu_t
  node_lxu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeLxuVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeLuuVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Luu_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Luu_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_luu();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeLuuVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Luu_t
  node_luu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeLuuVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeHVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::H_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::H_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_h();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeHVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::H_t
  node_h(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeHVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeHxVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Hx_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Hx_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_hx();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeHxVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Hx_t
  node_hx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeHxVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeHuVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Hu_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Hu_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_hu();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeHuVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Hu_t
  node_hu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeHuVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::G_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::G_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_g();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeGVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::G_t
  node_g(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeGVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGxVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Gx_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Gx_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_gx();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeGxVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Gx_t
  node_gx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeGxVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  struct NodeGuVisitor
      : boost::static_visitor<
            typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Gu_t>
  {
    using ReturnType = typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Gu_t;

    template <typename NodeDataDerived>
    ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
    {
      return node_data.node_gu();
    }

    static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
    {
      return boost::apply_visitor(NodeGuVisitor(), node_data);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
  inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Gu_t
  node_gu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
  {
    return NodeGuVisitor<Scalar, Options, NodeCollectionTpl>::run(node_data);
  }

} // namespace galileo

#endif // __galileo_core_node_basic_visitors_hxx__
