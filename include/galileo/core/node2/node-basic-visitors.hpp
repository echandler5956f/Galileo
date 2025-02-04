#ifndef __galileo_core_node_basic_visitors_hpp__
#define __galileo_core_node_basic_visitors_hpp__

#include "galileo/core/node/fwd.hpp"

namespace galileo
{

    // Visitors on Node Models

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline NodeDataTpl<Scalar, Options, NodeCollectionTpl>
    createData(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateVectorType, typename ControlVectorType>
    inline void calc(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                     NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data,
                     const Eigen::MatrixBase<StateVectorType> &xs,
                     const Eigen::MatrixBase<ControlVectorType> &us);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateVectorType, typename ControlVectorType>
    inline void calcDiff(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                         NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data,
                         const Eigen::MatrixBase<StateVectorType> &xs,
                         const Eigen::MatrixBase<ControlVectorType> &us);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::X_t
    node_get_x_lb(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::X_t
    node_get_x_ub(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::U_t
    node_get_u_lb(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::U_t
    node_get_u_ub(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::H_t
    node_get_h_lb(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::H_t
    node_get_h_ub(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::G_t
    node_get_g_lb(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeModelTpl<Scalar, Options, NodeCollectionTpl>::G_t
    node_get_g_ub(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline Eigen::Index
    node_get_nh(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline Eigen::Index
    node_get_ng(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline Eigen::Index
    node_get_id(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateBoundVectorType>
    inline void
    node_set_x_lb(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                  const Eigen::MatrixBase<StateBoundVectorType> &x_lb);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename StateBoundVectorType>
    inline void
    node_set_x_ub(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                  const Eigen::MatrixBase<StateBoundVectorType> &x_ub);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename ControlBoundVectorType>
    inline void
    node_set_u_lb(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                  const Eigen::MatrixBase<ControlBoundVectorType> &u_lb);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename ControlBoundVectorType>
    inline void
    node_set_u_ub(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                  const Eigen::MatrixBase<ControlBoundVectorType> &u_ub);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename EqualityBoundVectorType>
    inline void
    node_set_h_lb(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                  const Eigen::MatrixBase<EqualityBoundVectorType> &h_lb);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename EqualityBoundVectorType>
    inline void
    node_set_h_ub(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                  const Eigen::MatrixBase<EqualityBoundVectorType> &h_ub);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename InequalityBoundVectorType>
    inline void
    node_set_g_lb(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                  const Eigen::MatrixBase<InequalityBoundVectorType> &g_lb);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl, typename InequalityBoundVectorType>
    inline void
    node_set_g_ub(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                  const Eigen::MatrixBase<InequalityBoundVectorType> &g_ub);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline void
    node_set_nh(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                Eigen::Index nh);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline void
    node_set_ng(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                Eigen::Index ng);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline void
    node_set_id(NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                Eigen::Index id);

    template <typename NewScalar, typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    typename CastType<NewScalar, NodeModelTpl<Scalar, Options, NodeCollectionTpl>>::type
    cast_node(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model);

    // Visitors on Node Data

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::dX_t
    node_xdot(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Fx_t
    node_fx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Fu_t
    node_fu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::L_t
    node_l(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lx_t
    node_lx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lu_t
    node_lu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lxx_t
    node_lxx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Lxu_t
    node_lxu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Luu_t
    node_luu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::H_t
    node_h(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Hx_t
    node_hx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Hu_t
    node_hu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::G_t
    node_g(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Gx_t
    node_gx(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

    template <typename Scalar, int Options, template <typename S, int O> class NodeCollectionTpl>
    inline typename NodeDataTpl<Scalar, Options, NodeCollectionTpl>::Gu_t
    node_gu(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data);

} // namespace galileo

/* --- Details -------------------------------------------------------------------- */
// Included later
// #include "galileo/core/node/core-basic-visitors.hxx"

#endif // __galileo_core_node_basic_visitors_hpp__
