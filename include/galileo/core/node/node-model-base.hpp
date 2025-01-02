#ifndef __galileo_core_node_model_base_hpp__
#define __galileo_core_node_model_base_hpp__

#include "galileo/core/node/node-base.hpp"
#include <limits>

#define GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, TYPENAME)            \
    typedef TYPENAME traits<Node>::Scalar Scalar;                     \
    typedef TYPENAME traits<Node>::NodeModelDerived NodeModelDerived; \
    typedef TYPENAME traits<Node>::NodeDataDerived NodeDataDerived;   \
    enum                                                              \
    {                                                                 \
        Options = traits<Node>::Options,                              \
        NX = traits<Node>::NX,                                        \
        NU = traits<Node>::NU,                                        \
        NDX = traits<Node>::NDX,                                      \
        NH = traits<Node>::NH,                                        \
        NG = traits<Node>::NG                                         \
    };                                                                \
    typedef TYPENAME traits<Node>::X_t X_t;                           \
    typedef TYPENAME traits<Node>::U_t U_t;                           \
    typedef TYPENAME traits<Node>::dX_t dX_t;                         \
    typedef TYPENAME traits<Node>::Fx_t Fx_t;                         \
    typedef TYPENAME traits<Node>::Fu_t Fu_t;                         \
    typedef TYPENAME traits<Node>::L_t L_t;                           \
    typedef TYPENAME traits<Node>::Lx_t Lx_t;                         \
    typedef TYPENAME traits<Node>::Lu_t Lu_t;                         \
    typedef TYPENAME traits<Node>::Lxx_t Lxx_t;                       \
    typedef TYPENAME traits<Node>::Lxu_t Lxu_t;                       \
    typedef TYPENAME traits<Node>::Luu_t Luu_t;                       \
    typedef TYPENAME traits<Node>::H_t H_t;                           \
    typedef TYPENAME traits<Node>::Hx_t Hx_t;                         \
    typedef TYPENAME traits<Node>::Hu_t Hu_t;                         \
    typedef TYPENAME traits<Node>::G_t G_t;                           \
    typedef TYPENAME traits<Node>::Gx_t Gx_t;                         \
    typedef TYPENAME traits<Node>::Gu_t Gu_t;

#ifdef __clang__

#define GALILEO_NODE_TYPEDEF(Node) \
    GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, GALILEO_EMPTY_ARG)
#define GALILEO_NODE_TYPEDEF_TEMPLATE(Node) \
    GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, typename)

#elif (__GNUC__ == 4) && (__GNUC_MINOR__ == 4) && (__GNUC_PATCHLEVEL__ == 2)

#define GALILEO_NODE_TYPEDEF(Node) \
    GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, GALILEO_EMPTY_ARG)
#define GALILEO_NODE_TYPEDEF_TEMPLATE(Node) \
    GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, typename)

#else

#define GALILEO_NODE_TYPEDEF(Node) GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, typename)
#define GALILEO_NODE_TYPEDEF_TEMPLATE(Node) \
    GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, typename)

#endif

#define GALILEO_NODE_CAST_TYPE_SPECIALIZATION(NodeModelTpl) \
    template <typename Scalar, typename NewScalar>          \
    struct CastType<NewScalar, NodeModelTpl<Scalar>>        \
    {                                                       \
        typedef NodeModelTpl<NewScalar> type;               \
    }

namespace galileo
{
    template <typename Derived>
    class NodeModelBase : NumericalBase<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::NodeDerived NodeDerived;
        GALILEO_NODE_TYPEDEF_TEMPLATE(NodeDerived);

        NodeModelDerived &derived()
        {
            return *static_cast<Derived *>(this);
        }

        const NodeModelDerived &derived() const
        {
            return *static_cast<const Derived *>(this);
        }

        NodeDataDerived createData() const
        {
            return derived().createData();
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(NodeDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &xs,
                  const Eigen::MatrixBase<ControlVectorType> &us) const
        {
            derived().calc(data, xs.derived(), us.derived());
        }

        template <typename StateVectorType>
        void calc(NodeDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &xs) const
        {
            derived().calc(data, xs.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &xs,
                      const Eigen::MatrixBase<ControlVectorType> &us) const
        {
            derived().calcDiff(data, xs.derived(), us.derived());
        }

        template <typename StateVectorType>
        void calcDiff(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &xs) const
        {
            derived().calcDiff(data, xs.derived());
        }

        const X_t &get_x_lb() const
        {
            return derived().get_x_lb();
        }

        const X_t &get_x_ub() const
        {
            return derived().get_x_ub();
        }

        const U_t &get_u_lb() const
        {
            return derived().get_u_lb();
        }

        const U_t &get_u_ub() const
        {
            return derived().get_u_ub();
        }

        const H_t &get_h_lb() const
        {
            return derived().get_h_lb();
        }

        const H_t &get_h_ub() const
        {
            return derived().get_h_ub();
        }

        const G_t &get_g_lb() const
        {
            return derived().get_g_lb();
        }

        const G_t &get_g_ub() const
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
        void set_h_lb(const Eigen::MatrixBase<EqualityBoundVectorType> &h_lb)
        {
            derived().set_h_lb(h_lb);
        }

        template <typename EqualityBoundVectorType>
        void set_h_ub(const Eigen::MatrixBase<EqualityBoundVectorType> &h_ub)
        {
            derived().set_h_ub(h_ub);
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

        template <typename NewScalar>
        typename CastType<NewScalar, Derived>::type cast() const
        {
            return derived().template cast<NewScalar>();
        }

        // template <class OtherDerived>
        // bool operator==(const NodeModelBase<OtherDerived> &other) const
        // {
        //     return derived().isEqual(other.derived());
        // }

        // template <class OtherDerived>
        // bool operator!=(const NodeModelBase<OtherDerived> &other) const
        // {
        //     return !(internal::comparison_eq(derived(), other.derived()));
        // }

        // template <class OtherDerived>
        // bool isEqual(const NodeModelBase<OtherDerived> &) const
        // {
        //     return false;
        // }

        // bool isEqual(const NodeModelBase<Derived> &other) const
        // {
        //     return derived().isEqual(other.derived());
        // }

    protected:
        // Default constructor: protected.
        // Prevent the construction of stand-alone NodeModelBase.
        inline NodeModelBase() : nh_(std::numeric_limits<Eigen::Index>::max()), ng_(std::numeric_limits<Eigen::Index>::max(), id_(std::numeric_limits<Eigen::Index>::max()))
        {
        }

        // Copy constructor: protected.
        // Copy of stand-alone NodeModelBase are prevented, but can be used from inheriting
        // objects. Copy is done by calling copy operator.
        inline NodeModelBase(const NodeModelBase &clone)
        {
            *this = clone;
        }

        // Copy operator: protected.
        // Copy of stand-alone NodeModelBase are prevented, but can be used from inheriting
        // objects.
        inline NodeModelBase &operator=(const NodeModelBase &clone)
        {
            x_lb_ = clone.x_lb_;
            x_ub_ = clone.x_ub_;

            u_lb_ = clone.u_lb_;
            u_ub_ = clone.u_ub_;

            h_lb_ = clone.h_lb_;
            h_ub_ = clone.h_ub_;
            g_lb_ = clone.g_lb_;
            g_ub_ = clone.g_ub_;

            nh_ = clone.nh_;
            ng_ = clone.ng_;

            id_ = clone.id_;
            return *this;
        }

        X_t x_lb_; // Lower state limits
        X_t x_ub_; // Upper state limits

        U_t u_lb_; // Lower control limits
        U_t u_ub_; // Upper control limits

        H_t h_lb_; // Lower bound of the equality constraints
        H_t h_ub_; // Upper bound of the equality constraints
        G_t g_lb_; // Lower bound of the inequality constraints
        G_t g_ub_; // Upper bound of the inequality constraints

        Eigen::Index nh_; // Number of equality constraints
        Eigen::Index ng_; // Number of inequality constraints

        Eigen::Index id_; // Index of the node in the segment list

    }; // class NodeModelBase

} // namespace galileo

#endif // __galileo_core_node_model_base_hpp__
