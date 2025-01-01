#ifndef __galileo_core_node_data_base_hpp__
#define __galileo_core_node_data_base_hpp__

#include "galileo/core/node/node-base.hpp"
#include "galileo/core/node/node-model-base.hpp"

#define GALILEO_NODE_DATA_TYPEDEF_GENERIC(Node, TYPENAME)           \
    GALILEO_NODE_MODEL_TYPEDEF_GENERIC(Node, TYPENAME);             \
    typedef TYPENAME traits<Node>::XTypeConstRef XTypeConstRef;     \
    typedef TYPENAME traits<Node>::XTypeRef XTypeRef;               \
    typedef TYPENAME traits<Node>::dXTypeConstRef dXTypeConstRef;   \
    typedef TYPENAME traits<Node>::dXTypeRef dXTypeRef;             \
    typedef TYPENAME traits<Node>::UTypeConstRef UTypeConstRef;     \
    typedef TYPENAME traits<Node>::UTypeRef UTypeRef;               \
    typedef TYPENAME traits<Node>::FxTypeConstRef FxTypeConstRef;   \
    typedef TYPENAME traits<Node>::FxTypeRef FxTypeRef;             \
    typedef TYPENAME traits<Node>::FuTypeConstRef FuTypeConstRef;   \
    typedef TYPENAME traits<Node>::FuTypeRef FuTypeRef;             \
    typedef TYPENAME traits<Node>::LTypeConstRef LTypeConstRef;     \
    typedef TYPENAME traits<Node>::LTypeRef LTypeRef;               \
    typedef TYPENAME traits<Node>::LxTypeConstRef LxTypeConstRef;   \
    typedef TYPENAME traits<Node>::LxTypeRef LxTypeRef;             \
    typedef TYPENAME traits<Node>::LuTypeConstRef LuTypeConstRef;   \
    typedef TYPENAME traits<Node>::LuTypeRef LuTypeRef;             \
    typedef TYPENAME traits<Node>::LxxTypeConstRef LxxTypeConstRef; \
    typedef TYPENAME traits<Node>::LxxTypeRef LxxTypeRef;           \
    typedef TYPENAME traits<Node>::LxuTypeConstRef LxuTypeConstRef; \
    typedef TYPENAME traits<Node>::LxuTypeRef LxuTypeRef;           \
    typedef TYPENAME traits<Node>::LuuTypeConstRef LuuTypeConstRef; \
    typedef TYPENAME traits<Node>::LuuTypeRef LuuTypeRef;           \
    typedef TYPENAME traits<Node>::HTypeConstRef HTypeConstRef;     \
    typedef TYPENAME traits<Node>::HTypeRef HTypeRef;               \
    typedef TYPENAME traits<Node>::HxTypeConstRef HxTypeConstRef;   \
    typedef TYPENAME traits<Node>::HxTypeRef HxTypeRef;             \
    typedef TYPENAME traits<Node>::HuTypeConstRef HuTypeConstRef;   \
    typedef TYPENAME traits<Node>::HuTypeRef HuTypeRef;             \
    typedef TYPENAME traits<Node>::GTypeConstRef GTypeConstRef;     \
    typedef TYPENAME traits<Node>::GTypeRef GTypeRef;               \
    typedef TYPENAME traits<Node>::GxTypeConstRef GxTypeConstRef;   \
    typedef TYPENAME traits<Node>::GxTypeRef GxTypeRef;             \
    typedef TYPENAME traits<Node>::GuTypeConstRef GuTypeConstRef;   \
    typedef TYPENAME traits<Node>::GuTypeRef GuTypeRef;

#ifdef __clang__

#define GALILEO_NODE_DATA_TYPEDEF(Node) \
    GALILEO_NODE_DATA_TYPEDEF_GENERIC(Node, GALILEO_EMPTY_ARG)
#define GALILEO_NODE_DATA_TYPEDEF_TEMPLATE(Node) \
    GALILEO_NODE_DATA_TYPEDEF_GENERIC(Node, typename)

#elif (__GNUC__ == 4) && (__GNUC_MINOR__ == 4) && (__GNUC_PATCHLEVEL__ == 2)

#define GALILEO_NODE_DATA_TYPEDEF(Node) \
    GALILEO_NODE_DATA_TYPEDEF_GENERIC(Node, GALILEO_EMPTY_ARG)
#define GALILEO_NODE_DATA_TYPEDEF_TEMPLATE(Node) \
    GALILEO_NODE_DATA_TYPEDEF_GENERIC(Node, typename)

#else

#define GALILEO_NODE_DATA_TYPEDEF(Node) GALILEO_NODE_DATA_TYPEDEF_GENERIC(Node, typename)
#define GALILEO_NODE_DATA_TYPEDEF_TEMPLATE(Node) \
    GALILEO_NODE_DATA_TYPEDEF_GENERIC(Node, typename)

#endif

#define GALILEO_NODE_DATA_BASE_DEFAULT_ACCESSOR \
    dXTypeConstRef XDot_accessor() const        \
    {                                           \
        return XDot;                            \
    }                                           \
    dXTypeRef XDot_accessor()                   \
    {                                           \
        return XDot;                            \
    }                                           \
    FxTypeConstRef Fx_accessor() const          \
    {                                           \
        return Fx;                              \
    }                                           \
    FxTypeRef Fx_accessor()                     \
    {                                           \
        return Fx;                              \
    }                                           \
    FuTypeConstRef Fu_accessor() const          \
    {                                           \
        return Fu;                              \
    }                                           \
    FuTypeRef Fu_accessor()                     \
    {                                           \
        return Fu;                              \
    }                                           \
    LTypeConstRef L_accessor() const            \
    {                                           \
        return L;                               \
    }                                           \
    LTypeRef L_accessor()                       \
    {                                           \
        return L;                               \
    }                                           \
    LxTypeConstRef Lx_accessor() const          \
    {                                           \
        return Lx;                              \
    }                                           \
    LxTypetRef Lx_accessor()                    \
    {                                           \
        return Lx;                              \
    }                                           \
    LuTypeConstRef Lu_accessor() const          \
    {                                           \
        return Lu;                              \
    }                                           \
    LuTypeRef Lu_accessor()                     \
    {                                           \
        return Lu;                              \
    }                                           \
    LxxTypeConstRef Lxx_accessor() const        \
    {                                           \
        return Lxx;                             \
    }                                           \
    LxxTypeRef Lxx_accessor()                   \
    {                                           \
        return Lxx;                             \
    }                                           \
    LxuTypeConstRef Lxu_accessor() const        \
    {                                           \
        return Lxu;                             \
    }                                           \
    LxuTypeRef Lxu_accessor()                   \
    {                                           \
        return Lxu;                             \
    }                                           \
    LuuTypeConstRef Luu_accessor() const        \
    {                                           \
        return Luu;                             \
    }                                           \
    LuuTypeRef Luu_accessor()                   \
    {                                           \
        return Luu;                             \
    }                                           \
    HTypeConstRef H_accessor() const            \
    {                                           \
        return H;                               \
    }                                           \
    HTypeRef H_accessor()                       \
    {                                           \
        return H;                               \
    }                                           \
    HxTypeConstRef Hx_accessor() const          \
    {                                           \
        return Hx;                              \
    }                                           \
    HxTypeRef Hx_accessor()                     \
    {                                           \
        return Hx;                              \
    }                                           \
    HuTypeConstRef Hu_accessor() const          \
    {                                           \
        return Hu;                              \
    }                                           \
    HuTypeRef Hu_accessor()                     \
    {                                           \
        return Hu;                              \
    }                                           \
    GTypeConstRef G_accessor() const            \
    {                                           \
        return G;                               \
    }                                           \
    GTypeRef G_accessor()                       \
    {                                           \
        return G;                               \
    }                                           \
    GxTypeConstRef Gx_accessor() const          \
    {                                           \
        return Gx;                              \
    }                                           \
    GxTypeRef Gx_accessor()                     \
    {                                           \
        return Gx;                              \
    }                                           \
    GuTypeConstRef Gu_accessor() const          \
    {                                           \
        return Gu;                              \
    }                                           \
    GuTypeRef Gu_accessor()                     \
    {                                           \
        return Gu;                              \
    }

// dX_t XDot; // XDot at x, u
// MatrixXs_t Fx;   // Jacobian of the dynamics w.r.t. the state \f$\mathbf{x}\f$
// MatrixXs_t Fu;   // Jacobian of the dynamics w.r.t. the control \f$\mathbf{u}\f$

// Scalar cost;    // Cost at x, u
// VectorXs_t Lx;  // Jacobian of the cost w.r.t. the state \f$\mathbf{x}\f$
// VectorXs_t Lu;  // Jacobian of the cost w.r.t. the control \f$\mathbf{u}\f$
// MatrixXs_t Lxx; // Hessian of the cost w.r.t. the state \f$\mathbf{x}\f$
// MatrixXs_t Lxu; // Hessian of the cost w.r.t. the state \f$\mathbf{x}\f$ and
//                 // control u
// MatrixXs_t Luu; // Hessian of the cost w.r.t. the control \f$\mathbf{u}\f$

// VectorXs_t g;  // Inequality constraint values
// MatrixXs_t Gx; // Jacobian of the inequality constraint w.r.t. the state
//                // \f$\mathbf{x}\f$
// MatrixXs_t Gu; // Jacobian of the inequality constraint w.r.t. the control
//                // \f$\mathbf{u}\f$

// VectorXs_t h;  // Equality constraint values
// MatrixXs_t Hx; // Jacobian of the equality constraint w.r.t. the state
//                // \f$\mathbf{x}\f$
// MatrixXs_t Hu; // Jacobian of the equality constraint w.r.t the control
//                // \f$\mathbf{u}\f$

#define GALILEO_NODE_DATA_BASE_ACCESSOR_DEFAULT_RETURN_TYPE       \
    typedef const X_t &XTypeConstRef;                             \
    typedef X_t &XTypeRef;                                        \
    typedef const dX_t &dXTypeConstRef;                           \
    typedef dX_t &dXTypeRef;                                      \
    typedef const U_t &UTypeConstRef;                             \
    typedef U_t &UTypeRef;                                        \
    tyedef const H_t &HTypeConstRef;                              \
    typedef H_t &HTypeRef;                                        \
    typedef const G_t &GTypeConstRef;                             \
    typedef G_t &GTypeRef;                                        \
    typedef const Eigen::Matrix<Scalar, NDX, NX> &FxTypeConstRef; \
    typedef Eigen::Matrix<Scalar, NDX, NX> &FxTypeRef;            \
    typedef const Eigen::Matrix<Scalar, NDX, NU> &FuTypeConstRef; \
    typedef Eigen::Matrix<Scalar, NDX, NU> &FuTypeRef;            \
    typedef const Scalar &LTypeConstRef;                          \
    typedef Scalar &LTypeRef;                                     \
    typedef const Eigen::Matrix<Scalar, 1, NX> &LxTypeConstRef;   \
    typedef Eigen::Matrix<Scalar, 1, NX> &LxTypeRef;              \
    typedef const Eigen::Matrix<Scalar, 1, NU> &LuTypeConstRef;   \
    typedef Eigen::Matrix<Scalar, 1, NU> &LuTypeRef;              \
    typedef const Eigen::Matrix<Scalar, NX, NX> &LxxTypeConstRef; \
    typedef Eigen::Matrix<Scalar, NX, NX> &LxxTypeRef;            \
    typedef const Eigen::Matrix<Scalar, NX, NU> &LxuTypeConstRef; \
    typedef Eigen::Matrix<Scalar, NX, NU> &LxuTypeRef;            \
    typedef const Eigen::Matrix<Scalar, NU, NU> &LuuTypeConstRef; \
    typedef Eigen::Matrix<Scalar, NU, NU> &LuuTypeRef;            \
    typedef const Eigen::Matrix<Scalar, NH, NX> &HxTypeConstRef;  \
    typedef Eigen::Matrix<Scalar, NH, NX> &HxTypeRef;             \
    typedef const Eigen::Matrix<Scalar, NH, NU> &HuTypeConstRef;  \
    typedef Eigen::Matrix<Scalar, NH, NU> &HuTypeRef;             \
    typedef const Eigen::Matrix<Scalar, NG, NX> &GxTypeConstRef;  \
    typedef Eigen::Matrix<Scalar, NG, NX> &GxTypeRef;             \
    typedef const Eigen::Matrix<Scalar, NG, NU> &GuTypeConstRef;  \
    typedef Eigen::Matrix<Scalar, NG, NU> &GuTypeRef;

namespace galileo
{

    template <typename Derived>
    struct NodeDataBase : NumericalBase<Derived>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::NodeDerived NodeDerived;
        GALILEO_NODE_DATA_TYPEDEF_TEMPLATE(NodeDerived);

        Derived &derived()
        {
            return *static_cast<Derived *>(this);
        }
        const Derived &derived() const
        {
            return *static_cast<const Derived *>(this);
        }

        dXTypeConstRef XDot() const
        {
            return derived().XDot_accessor();
        }
        dXTypeRef XDot()
        {
            return derived().XDot_accessor();
        }

        FxTypeConstRef Fx() const
        {
            return derived().Fx_accessor();
        }
        FxTypeRef Fx()
        {
            return derived().Fx_accessor();
        }

        FuTypeConstRef Fu() const
        {
            return derived().Fu_accessor();
        }
        FuTypeRef Fu()
        {
            return derived().Fu_accessor();
        }

        LTypeConstRef L() const
        {
            return derived().L_accessor();
        }
        LTypeRef L()
        {
            return derived().L_accessor();
        }

        LxTypeConstRef Lx() const
        {
            return derived().Lx_accessor();
        }
        LxTypeRef Lx()
        {
            return derived().Lx_accessor();
        }

        LuTypeConstRef Lu() const
        {
            return derived().Lu_accessor();
        }
        LuTypeRef Lu()
        {
            return derived().Lu_accessor();
        }

        LxxTypeConstRef Lxx() const
        {
            return derived().Lxx_accessor();
        }
        LxxTypeRef Lxx()
        {
            return derived().Lxx_accessor();
        }

        LxuTypeConstRef Lxu() const
        {
            return derived().Lxu_accessor();
        }
        LxuTypeRef Lxu()
        {
            return derived().Lxu_accessor();
        }

        LuuTypeConstRef Luu() const
        {
            return derived().Luu_accessor();
        }
        LuuTypeRef Luu()
        {
            return derived().Luu_accessor();
        }

        HTypeConstRef H() const
        {
            return derived().H_accessor();
        }
        HTypeRef H()
        {
            return derived().H_accessor();
        }

        HxTypeConstRef Hx() const
        {
            return derived().Hx_accessor();
        }
        HxTypeRef Hx()
        {
            return derived().Hx_accessor();
        }

        HuTypeConstRef Hu() const
        {
            return derived().Hu_accessor();
        }
        HuTypeRef Hu()
        {
            return derived().Hu_accessor();
        }

        GTypeConstRef G() const
        {
            return derived().G_accessor();
        }
        GTypeRef G()
        {
            return derived().G_accessor();
        }

        GxTypeConstRef Gx() const
        {
            return derived().Gx_accessor();
        }
        GxTypeRef Gx()
        {
            return derived().Gx_accessor();
        }

        GuTypeConstRef Gu() const
        {
            return derived().Gu_accessor();
        }
        GuTypeRef Gu()
        {
            return derived().Gu_accessor();
        }

        // template <typename OtherDerived>
        // bool operator==(const NodeDataBase<OtherDerived> &other) const
        // {
        //     return derived().isEqual(other.derived());
        // }

        // // Default operator== implementation
        // bool isEqual(const NodeDataBase<Derived> & other) const
        // {
        // return internal::comparison_eq(joint_q(), other.joint_q())
        //         && internal::comparison_eq(joint_v(), other.joint_v())
        //         && internal::comparison_eq(S(), other.S()) && internal::comparison_eq(M(), other.M())
        //         && internal::comparison_eq(v(), other.v()) && internal::comparison_eq(c(), other.c())
        //         && internal::comparison_eq(U(), other.U())
        //         && internal::comparison_eq(Dinv(), other.Dinv())
        //         && internal::comparison_eq(UDinv(), other.UDinv());
        // }

        // // Default operator== implementation
        // template <typename OtherDerived>
        // bool isEqual(const NodeDataBase<OtherDerived> & /*other*/) const
        // {
        //     return false;
        // }

        // bool operator!=(const NodeDataBase<Derived> &other) const
        // {
        //     return derived().isNotEqual(other.derived());
        // }

        // // Default operator!= implementation
        // bool isNotEqual(const NodeDataBase<Derived> &other) const
        // {
        //     return !(internal::comparison_eq(derived(), other.derived()));
        // }

    protected:
        // Default constructor: protected.
        inline NodeDataBase()
        {
        }

    }; // struct NodeDataBase

} // namespace galileo

#endif // __galileo_core_node_data_base_hpp__