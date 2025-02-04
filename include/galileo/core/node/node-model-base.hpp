#ifndef __galileo_core_node_model_base_hpp__
#define __galileo_core_node_model_base_hpp__

#include "galileo/core/node/node-base.hpp"
#include <memory>

#define GALILEO_NODE_BASIC_TYPEDEF(Node)                              \
    using Scalar = typename traits<Node>::Scalar;                     \
    using NodeModelDerived = typename traits<Node>::NodeModelDerived; \
    using NodeDataDerived = typename traits<Node>::NodeDataDerived;

#define GALILEO_NODE_ACTION_DEF_TYPEDEF(Node)                                     \
    using State_t = typename traits<Node>::State_t;                               \
    using Actuation_t = typename traits<Node>::Actuation_t;                       \
    using ConstraintModelCollection_t = typename traits<Node>::ConstraintModelCollection_t; \
    using CostModelCollection_t = typename traits<Node>::CostModelCollection_t;

#define GALILEO_NODE_MODEL_TYPEDEF(Node)                  \
    using X_Bounds_t = typename traits<Node>::X_Bounds_t; \
    using U_Bounds_t = typename traits<Node>::U_Bounds_t; \
    using H_Bounds_t = typename traits<Node>::H_Bounds_t; \
    using G_Bounds_t = typename traits<Node>::G_Bounds_t;

#define GALILEO_NODE_DATA_TYPEDEF(Node)         \
    using F_t = typename traits<Node>::F_t;     \
    using Fx_t = typename traits<Node>::Fx_t;   \
    using Fu_t = typename traits<Node>::Fu_t;   \
    using L_t = typename traits<Node>::L_t;     \
    using Lx_t = typename traits<Node>::Lx_t;   \
    using Lu_t = typename traits<Node>::Lu_t;   \
    using Lxx_t = typename traits<Node>::Lxx_t; \
    using Lxu_t = typename traits<Node>::Lxu_t; \
    using Luu_t = typename traits<Node>::Luu_t; \
    using H_t = typename traits<Node>::H_t;     \
    using Hx_t = typename traits<Node>::Hx_t;   \
    using Hu_t = typename traits<Node>::Hu_t;   \
    using G_t = typename traits<Node>::G_t;     \
    using Gx_t = typename traits<Node>::Gx_t;   \
    using Gu_t = typename traits<Node>::Gu_t;

namespace galileo
{

    template <typename Derived>
    class NodeModelBase : CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using NodeDerived = typename traits<Derived>::NodeDerived;
        GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
        GALILEO_NODE_ACTION_DEF_TYPEDEF(NodeDerived);
        GALILEO_NODE_MODEL_TYPEDEF(NodeDerived);

        NodeDataDerived createData() const
        {
            return derived().createData();
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(NodeDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            derived().calc(data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calc(NodeDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            derived().calc(data, x.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            derived().calcDiff(data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calcDiff(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            derived().calcDiff(data, x.derived());
        }

        const X_Bounds_t &get_x_lb() const
        {
            return derived().get_x_lb();
        }

        const X_Bounds_t &get_x_ub() const
        {
            return derived().get_x_ub();
        }

        const U_Bounds_t &get_u_lb() const
        {
            return derived().get_u_lb();
        }

        const U_Bounds_t &get_u_ub() const
        {
            return derived().get_u_ub();
        }

        const H_Bounds_t &get_h_eq() const
        {
            return derived().get_h_eq();
        }

        const G_Bounds_t &get_g_lb() const
        {
            return derived().get_g_lb();
        }

        const G_Bounds_t &get_g_ub() const
        {
            return derived().get_g_ub();
        }

        Eigen::Index get_nh() const
        {
            return derived().get_nh();
        }

        Eigen::Index get_ng() const
        {
            return derived().get_ng();
        }

        Eigen::Index get_id() const
        {
            return derived().get_id();
        }

        template <typename StateBoundVectorType>
        void set_x_lb(const Eigen::MatrixBase<StateBoundVectorType> &x_lb)
        {
            derived().set_x_lb(x_lb);
        }

        template <typename StateBoundVectorType>
        void set_x_ub(const Eigen::MatrixBase<StateBoundVectorType> &x_ub)
        {
            derived().set_x_ub(x_ub);
        }

        template <typename ControlBoundVectorType>
        void set_u_lb(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb)
        {
            derived().set_u_lb(u_lb);
        }

        template <typename ControlBoundVectorType>
        void set_u_ub(const Eigen::MatrixBase<ControlBoundVectorType> &u_ub)
        {
            derived().set_u_ub(u_ub);
        }

        template <typename EqualityBoundVectorType>
        void set_h_eq(const Eigen::MatrixBase<EqualityBoundVectorType> &h_eq)
        {
            derived().set_h_eq(h_eq);
        }

        template <typename InequalityBoundVectorType>
        void set_g_lb(const Eigen::MatrixBase<InequalityBoundVectorType> &g_lb)
        {
            derived().set_g_lb(g_lb);
        }

        template <typename InequalityBoundVectorType>
        void set_g_ub(const Eigen::MatrixBase<InequalityBoundVectorType> &g_ub)
        {
            derived().set_g_ub(g_ub);
        }

        void set_nh(Eigen::Index nh)
        {
            derived().set_nh(nh);
        }

        void set_ng(Eigen::Index ng)
        {
            derived().set_ng(ng);
        }

        void set_id(Eigen::Index id)
        {
            derived().set_id(id);
        }

    protected:
        inline NodeModelBase()
        {
        }

        inline NodeModelBase(const NodeModelBase &clone)
        {
            *this = clone;
        }

        inline NodeModelBase &operator=(const NodeModelBase &clone)
        {
            x_lb_ = clone.x_lb_;
            x_ub_ = clone.x_ub_;

            u_lb_ = clone.u_lb_;
            u_ub_ = clone.u_ub_;

            h_eq_ = clone.h_eq_;
            g_lb_ = clone.g_lb_;
            g_ub_ = clone.g_ub_;

            nh_ = clone.nh_;
            ng_ = clone.ng_;

            id_ = clone.id_;
            return *this;
        }

        X_Bounds_t x_lb_; // Lower state limits
        X_Bounds_t x_ub_; // Upper state limits

        U_Bounds_t u_lb_; // Lower control limits
        U_Bounds_t u_ub_; // Upper control limits

        H_Bounds_t h_eq_; // Lower bound of the equality constraints
        G_Bounds_t g_lb_; // Lower bound of the inequality constraints
        G_Bounds_t g_ub_; // Upper bound of the inequality constraints

        Eigen::Index nh_; // Number of equality constraints
        Eigen::Index ng_; // Number of inequality constraints

        Eigen::Index id_; // Index of the node in the segment list

    }; // class NodeModelBase
    
} // namespace galileo

#endif // __galileo_core_node_model_base_hpp__
