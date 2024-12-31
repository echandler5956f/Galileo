#ifndef __galileo_core_node_model_base_hpp__
#define __galileo_core_node_model_base_hpp__

#include "galileo/core/node/node-base.hpp"
#include <limits>

#define GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, TYPENAME)            \
    typedef TYPENAME traits<Node>::Scalar Scalar;                     \
    typedef TYPENAME traits<Node>::NodeModelDerived NodeModelDerived; \
    typedef TYPENAME traits<Node>::NodeDataDerived NodeDataDerived;   \
    typedef TYPENAME traits<Node>::State_t State_t;                   \
    typedef TYPENAME traits<Node>::Control_t Control_t;               \
    enum                                                              \
    {                                                                 \
        Options = traits<Node>::Options,                              \
        NX = traits<Node>::NX,                                        \
        NDX = traits<Node>::NDX,                                      \
        NU = traits<Node>::NU,                                        \
        NH = traits<Node>::NH,                                        \
        NG = traits<Node>::NG                                         \
    };                                                                \
    typedef TYPENAME traits<Node>::X_t X_t;                           \
    typedef TYPENAME traits<Node>::dX_t dX_t;                         \
    typedef TYPENAME traits<Node>::U_t U_t;                           \
    typedef TYPENAME traits<Node>::H_t H_t;                           \
    typedef TYPENAME traits<Node>::G_t G_t;

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

        template <typename StateVectorType, typename ControlVectorType>
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

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &xs) const
        {
            derived().calcDiff(data, xs.derived());
        }

        X_t &get_x_lb()
        {
            return derived().get_x_lb();
        }

        X_t &get_x_ub()
        {
            return derived().get_x_ub();
        }

        U_t &get_u_lb()
        {
            return derived().get_u_lb();
        }

        U_t &get_u_ub()
        {
            return derived().get_u_ub();
        }

        H_t &get_h_lb()
        {
            return derived().get_h_lb();
        }

        H_t &get_h_ub()
        {
            return derived().get_h_ub();
        }

        G_t &get_g_lb()
        {
            return derived().get_g_lb();
        }

        G_t &get_g_ub()
        {
            return derived().get_g_ub();
        }

        template <typename StateVectorType>
        void set_x_lb(const Eigen::MatrixBase<StateVectorType> &x_lb)
        {
            derived().set_x_lb(x_lb);
        }

        template <typename StateVectorType>
        void set_x_ub(const Eigen::MatrixBase<StateVectorType> &x_ub)
        {
            derived().set_x_ub(x_ub);
        }

        template <typename ControlVectorType>
        void set_u_lb(const Eigen::MatrixBase<ControlVectorType> &u_lb)
        {
            derived().set_u_lb(u_lb);
        }

        template <typename ControlVectorType>
        void set_u_ub(const Eigen::MatrixBase<ControlVectorType> &u_ub)
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
        inline NodeModelBase() : nh_(std::numeric_limits<std::size_t>::max()), ng_(std::numeric_limits<std::size_t>::max())
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

        std::size_t nh_; // Number of equality constraints
        std::size_t ng_; // Number of inequality constraints

    }; // class NodeModelBase

} // namespace galileo

#endif // __galileo_core_node_model_base_hpp__
