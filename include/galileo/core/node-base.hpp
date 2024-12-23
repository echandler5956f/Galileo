#pragma once

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <stdexcept>

#include "galileo/fwd.hpp"
#include "galileo/core/state-base.hpp"
#include "galileo/core/control-base.hpp"

#define GALILEO_NODE_TYPEDEF_GENERIC(Node, TYPENAME)                  \
    typedef TYPENAME traits<Node>::Scalar Scalar;                     \
    typedef TYPENAME traits<Node>::NodeModelDerived NodeModelDerived; \
    typedef TYPENAME traits<Node>::NodeDataDerived NodeDataDerived;   \
    typedef TYPENAME traits<Node>::StateModel StateModel;             \
    typedef TYPENAME traits<Node>::MathBase MathBase;                 \
    typedef TYPENAME MathBase::VectorXs VectorXs;                     \
    typedef TYPENAME MathBase::MatrixXs MatrixXs;

#define GALILEO_NODE_TYPEDEF_TEMPLATE(Node) \
    GALILEO_NODE_TYPEDEF_GENERIC(Node, typename)

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

        void calc(const boost::shared_ptr<NodeDataDerived> &data,
                  const Eigen::Ref<const VectorXs> &x,
                  const Eigen::Ref<const VectorXs> &u)
        {
            derived().calc(data, x, u);
        }

        void calc(const boost::shared_ptr<NodeDataDerived> &data,
                  const Eigen::Ref<const VectorXs> &x)
        {
            derived().calc(data, x);
        }

        void calcDiff(const boost::shared_ptr<NodeDataDerived> &data,
                      const Eigen::Ref<const VectorXs> &x,
                      const Eigen::Ref<const VectorXs> &u)
        {
            derived().calcDiff(data, x, u);
        }

        void calcDiff(const boost::shared_ptr<NodeDataDerived> &data,
                      const Eigen::Ref<const VectorXs> &x)
        {
            derived().calcDiff(data, x);
        }

        boost::shared_ptr<NodeDataDerived> createData() const
        {
            return derived().createData();
        }

        bool checkData(const boost::shared_ptr<NodeDataDerived> &data) const
        {
            return derived().checkData(data);
        }

        void quasiStatic(const boost::shared_ptr<NodeDataDerived> &data,
                         Eigen::Ref<VectorXs> u, const Eigen::Ref<const VectorXs> &x,
                         const std::size_t maxiter = 100, const Scalar tol = Scalar(1e-9))
        {
            derived().quasiStatic(data, u, x, maxiter, tol);
        }

        VectorXs quasiStatic_x(const boost::shared_ptr<NodeDataDerived> &data,
                               const VectorXs &x, const std::size_t maxiter = 100,
                               const Scalar tol = Scalar(1e-9))
        {
            derived().quasiStatic(data, x, maxiter, tol);
        }

        std::size_t get_nu() const
        {
            return derived().get_nu();
        }

        std::size_t get_nr() const
        {
            return derived().get_nr();
        }

        std::size_t get_ng() const
        {
            return derived().get_ng();
        }

        std::size_t get_nh() const
        {
            return derived().get_nh();
        }

        const boost::shared_ptr<StateModel> &get_state() const
        {
            return derived().get_state();
        }

        const VectorXs &get_g_lb() const
        {
            return derived().get_g_lb();
        }

        const VectorXs &get_g_ub() const
        {
            return derived().get_g_ub();
        }

        const VectorXs &get_u_lb() const
        {
            return derived().get_u_lb();
        }

        const VectorXs &get_u_ub() const
        {
            return derived().get_u_ub();
        }

        bool get_has_control_limits() const
        {
            return derived().get_has_control_limits();
        }

        void set_g_lb(const VectorXs &g_lb)
        {
            derived().set_g_lb(g_lb);
        }

        void set_g_ub(const VectorXs &g_ub)
        {
            derived().set_g_ub(g_ub);
        }

        void set_u_lb(const VectorXs &u_lb)
        {
            derived().set_u_lb(u_lb);
        }

        void set_u_ub(const VectorXs &u_ub)
        {
            derived().set_u_ub(u_ub);
        }

    protected:
        void update_has_control_limits()
        {
            derived().update_has_control_limits();
        }

        // Default constructor: protected.
        // Prevent the construction of stand-alone NodeModelBase.
        inline NodeModelBase() : nu_(0), nr_(0), ng_(0), nh_(0), has_control_limits_(false)
        {
            unone.resize(0);
            g_lb_.resize(0);
            g_ub_.resize(0);
            u_lb_.resize(0);
            u_ub_.resize(0);
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
            nu_ = clone.nu_;
            nr_ = clone.nr_;
            ng_ = clone.ng_;
            nh_ = clone.nh_;
            state_ = clone.state_;
            unone_ = clone.unone_;
            g_lb_ = clone.g_lb_;
            g_ub_ = clone.g_ub_;
            u_lb_ = clone.u_lb_;
            u_ub_ = clone.u_ub_;
            has_control_limits_ = clone.has_control_limits_;
            return *this;
        }

        std::size_t nu_;                      //!< Control dimension
        std::size_t nr_;                      //!< Dimension of the cost residual
        std::size_t ng_;                      //!< Number of inequality constraints
        std::size_t nh_;                      //!< Number of equality constraints
        boost::shared_ptr<StateModel> state_; //!< Model of the state
        VectorXs unone_;                      //!< Neutral state
        VectorXs g_lb_;                       //!< Lower bound of the inequality constraints
        VectorXs g_ub_;                       //!< Lower bound of the inequality constraints
        VectorXs u_lb_;                       //!< Lower control limits
        VectorXs u_ub_;                       //!< Upper control limits
        bool has_control_limits_;             //!< Indicates whether any of the control limits is
                                              //!< finite
    };

    template <typename Derived>
    struct NodeDataBase : NumericalBase<Derived>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::NodeDerived NodeDerived;
        GALILEO_NODE_TYPEDEF_TEMPLATE(NodeDerived);

        VectorXs xout; //!< xdot at x, u
        MatrixXs Fx;   //!< Jacobian of the dynamics w.r.t. the state \f$\mathbf{x}\f$
        MatrixXs Fu;   //!< Jacobian of the dynamics w.r.t. the control \f$\mathbf{u}\f$

        Scalar cost;  //!< Cost at x, u
        VectorXs r;   //!< Cost residual
        VectorXs Lx;  //!< Jacobian of the cost w.r.t. the state \f$\mathbf{x}\f$
        VectorXs Lu;  //!< Jacobian of the cost w.r.t. the control \f$\mathbf{u}\f$
        MatrixXs Lxx; //!< Hessian of the cost w.r.t. the state \f$\mathbf{x}\f$
        MatrixXs Lxu; //!< Hessian of the cost w.r.t. the state \f$\mathbf{x}\f$ and
                      //!< control u
        MatrixXs Luu; //!< Hessian of the cost w.r.t. the control \f$\mathbf{u}\f$

        VectorXs g;  //!< Inequality constraint values
        MatrixXs Gx; //!< Jacobian of the inequality constraint w.r.t. the state
                     //!< \f$\mathbf{x}\f$
        MatrixXs Gu; //!< Jacobian of the inequality constraint w.r.t. the control
                     //!< \f$\mathbf{u}\f$

        VectorXs h;  //!< Equality constraint values
        MatrixXs Hx; //!< Jacobian of the equality constraint w.r.t. the state
                     //!< \f$\mathbf{x}\f$
        MatrixXs Hu; //!< Jacobian of the equality constraint w.r.t the control
                     //!< \f$\mathbf{u}\f$

    protected:
        // Default constructor: protected.
        inline NodeDataBase()
        {
        }
    };

}