#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "galileo/fwd.hpp"

#define GALILEO_STATE_MODEL_TYPEDEF_GENERIC(State, TYPENAME)             \
    typedef TYPENAME traits<State>::Scalar Scalar;                       \
    typedef TYPENAME traits<State>::StateModelDerived StateModelDerived; \
    typedef TYPENAME traits<State>::VectorXs_t VectorXs_t;               \
    typedef TYPENAME traits<State>::MatrixXs_t MatrixXs_t;

#define GALILEO_STATE_TYPEDEF(State) GALILEO_STATE_MODEL_TYPEDEF_GENERIC(State, typename)
#define GALILEO_STATE_TYPEDEF_TEMPLATE(State) \
    GALILEO_STATE_MODEL_TYPEDEF_GENERIC(State, typename)

#define GALILEO_STATE_CAST_TYPE_SPECIALIZATION(StateModelTpl) \
    template <typename Scalar, typename NewScalar>            \
    struct CastType<NewScalar, StateModelTpl<Scalar>>         \
    {                                                         \
        typedef StateModelTpl<NewScalar> type;                \
    }

namespace galileo
{
    enum Jcomponent
    {
        both = 0,
        first = 1,
        second = 2
    };

    inline bool is_a_Jcomponent(Jcomponent firstsecond)
    {
        return (firstsecond == first || firstsecond == second || firstsecond == both);
    }

    /**
     * @brief Abstract class for the state representation
     *
     * A state is represented by its operators: difference, integrates, transport
     * and their derivatives. The difference operator returns the value of
     * \f$\mathbf{x}_{1}\ominus\mathbf{x}_{0}\f$ operation. Instead the integrate
     * operator returns the value of \f$\mathbf{x}\oplus\delta\mathbf{x}\f$. These
     * operators are used to compared two points on the state manifold
     * \f$\mathcal{M}\f$ or to advance the state given a tangential velocity
     * (\f$T_\mathbf{x} \mathcal{M}\f$). Therefore the points \f$\mathbf{x}\f$,
     * \f$\mathbf{x}_{0}\f$ and \f$\mathbf{x}_{1}\f$ belong to the manifold
     * \f$\mathcal{M}\f$; and \f$\delta\mathbf{x}\f$ or
     * \f$\mathbf{x}_{1}\ominus\mathbf{x}_{0}\f$ lie on its tangential space.
     */
    template <typename Derived>
    class StateBase : CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::StateDerived StateDerived;
        GALILEO_STATE_TYPEDEF_TEMPLATE(StateDerived);

        StateDerived &derived()
        {
            return *static_cast<Derived *>(this);
        }

        const StateDerived &derived() const
        {
            return *static_cast<const Derived *>(this);
        }

        /**
         * @brief Generate a zero state
         */
        VectorXs_t zero() const
        {
            derived().zero();
        }

        /**
         * @brief Generate a random state
         */
        VectorXs_t rand() const
        {
            derived().rand();
        }

        /**
         * @brief Compute the state manifold differentiation.
         *
         * The state differentiation is defined as:
         * \f{equation*}{
         *   \delta\mathbf{x} = \mathbf{x}_{1} \ominus \mathbf{x}_{0},
         * \f}
         * where \f$\mathbf{x}_{1}\f$, \f$\mathbf{x}_{0}\f$ are the current and
         * previous state which lie in a manifold \f$\mathcal{M}\f$, and
         * \f$\delta\mathbf{x} \in T_\mathbf{x} \mathcal{M}\f$ is the rate of change
         * in the state in the tangent space of the manifold.
         *
         * @param[in]  x0     Previous state point (size `nx`)
         * @param[in]  x1     Current state point (size `nx`)
         * @param[out] dxout  Difference between the current and previous state points
         * (size `ndx`)
         */
        void diff(const Eigen::Ref<const VectorXs_t> &x0,
                  const Eigen::Ref<const VectorXs_t> &x1,
                  Eigen::Ref<VectorXs_t> dxout) const
        {
            derived().diff(x0, x1, dxout);
        }

        /**
         * @brief Compute the state manifold integration.
         *
         * The state integration is defined as:
         * \f{equation*}{
         *   \mathbf{x}_{next} = \mathbf{x} \oplus \delta\mathbf{x},
         * \f}
         * where \f$\mathbf{x}\f$, \f$\mathbf{x}_{next}\f$ are the current and next
         * state which lie in a manifold \f$\mathcal{M}\f$, and \f$\delta\mathbf{x}
         * \in T_\mathbf{x} \mathcal{M}\f$ is the rate of change in the state in the
         * tangent space of the manifold.
         *
         * @param[in]  x     State point (size `nx`)
         * @param[in]  dx    Velocity vector (size `ndx`)
         * @param[out] xout  Next state point (size `nx`)
         */
        void integrate(const Eigen::Ref<const VectorXs_t> &x,
                       const Eigen::Ref<const VectorXs_t> &dx,
                       Eigen::Ref<VectorXs_t> xout) const
        {
            derived().integrate(x, dx, xout);
        }

        /**
         * @brief Compute the Jacobian of the state manifold differentiation.
         *
         * The state differentiation is defined as:
         * \f{equation*}{
         *   \delta\mathbf{x} = \mathbf{x}_{1} \ominus \mathbf{x}_{0},
         * \f}
         * where \f$\mathbf{x}_{1}\f$, \f$\mathbf{x}_{0}\f$ are the current and
         * previous state which lie in a manifold \f$\mathcal{M}\f$, and
         * \f$\delta\mathbf{x} \in T_\mathbf{x} \mathcal{M}\f$ is the rate of change
         * in the state in the tangent space of the manifold.
         *
         * The Jacobians lie in the tangent space of manifold, i.e.
         * \f$\mathbb{R}^{\textrm{ndx}\times\textrm{ndx}}\f$. Note that the state is
         * represented as a tuple of `nx` values and its dimension is `ndx`. Calling
         * \f$\boldsymbol{\Delta}(\mathbf{x}_{0}, \mathbf{x}_{1}) \f$, the difference
         * function, these Jacobians satisfy the following relationships:
         *  -
         * \f$\boldsymbol{\Delta}(\mathbf{x}_{0},\mathbf{x}_{0}\oplus\delta\mathbf{y})
         * - \boldsymbol{\Delta}(\mathbf{x}_{0},\mathbf{x}_{1}) =
         * \mathbf{J}_{\mathbf{x}_{1}}\delta\mathbf{y} +
         * \mathbf{o}(\mathbf{x}_{0})\f$.
         *  -
         * \f$\boldsymbol{\Delta}(\mathbf{x}_{0}\oplus\delta\mathbf{y},\mathbf{x}_{1})
         * - \boldsymbol{\Delta}(\mathbf{x}_{0},\mathbf{x}_{1}) =
         * \mathbf{J}_{\mathbf{x}_{0}}\delta\mathbf{y} +
         * \mathbf{o}(\mathbf{x}_{0})\f$,
         *
         * where \f$\mathbf{J}_{\mathbf{x}_{1}}\f$ and
         * \f$\mathbf{J}_{\mathbf{x}_{0}}\f$ are the Jacobian with respect to the
         * current and previous state, respectively.
         *
         * @param[in] x0           Previous state point (size `nx`)
         * @param[in] x1           Current state point (size `nx`)
         * @param[out] Jfirst      Jacobian of the difference operation relative to
         * the previous state point (size `ndx`\f$\times\f$`ndx`)
         * @param[out] Jsecond     Jacobian of the difference operation relative to
         * the current state point (size `ndx`\f$\times\f$`ndx`)
         * @param[in] firstsecond  Argument (either x0 and / or x1) with respect to
         * which the differentiation is performed.
         */
        void Jdiff(const Eigen::Ref<const VectorXs_t> &x0,
                   const Eigen::Ref<const VectorXs_t> &x1,
                   Eigen::Ref<MatrixXs_t> Jfirst, Eigen::Ref<MatrixXs_t> Jsecond,
                   const Jcomponent firstsecond = both) const
        {
            derived().Jdiff(x0, x1, Jfirst, Jsecond, firstsecond);
        }

        /**
         * @brief Compute the Jacobian of the state manifold integration.
         *
         * The state integration is defined as:
         * \f{equation*}{
         *   \mathbf{x}_{next} = \mathbf{x} \oplus \delta\mathbf{x},
         * \f}
         * where \f$\mathbf{x}\f$, \f$\mathbf{x}_{next}\f$ are the current and next
         * state which lie in a manifold \f$\mathcal{M}\f$, and \f$\delta\mathbf{x}
         * \in T_\mathbf{x} \mathcal{M}\f$ is the rate of change in the state in the
         * tangent space of the manifold.
         *
         * The Jacobians lie in the tangent space of manifold, i.e.
         * \f$\mathbb{R}^{\textrm{ndx}\times\textrm{ndx}}\f$. Note that the state is
         * represented as a tuple of `nx` values and its dimension is `ndx`. Calling
         * \f$ \mathbf{f}(\mathbf{x}, \delta\mathbf{x}) \f$, the integrate function,
         * these Jacobians satisfy the following relationships:
         *  -
         * \f$\mathbf{f}(\mathbf{x}\oplus\delta\mathbf{y},\delta\mathbf{x})\ominus\mathbf{f}(\mathbf{x},\delta\mathbf{x})
         * = \mathbf{J}_\mathbf{x}\delta\mathbf{y} + \mathbf{o}(\delta\mathbf{x})\f$.
         *  -
         * \f$\mathbf{f}(\mathbf{x},\delta\mathbf{x}+\delta\mathbf{y})\ominus\mathbf{f}(\mathbf{x},\delta\mathbf{x})
         * = \mathbf{J}_{\delta\mathbf{x}}\delta\mathbf{y} +
         * \mathbf{o}(\delta\mathbf{x})\f$,
         *
         * where \f$\mathbf{J}_{\delta\mathbf{x}}\f$ and \f$\mathbf{J}_{\mathbf{x}}\f$
         * are the Jacobian with respect to the state and velocity, respectively.
         *
         * @param[in] x            State point (size `nx`)
         * @param[in] dx           Velocity vector (size `ndx`)
         * @param[out] Jfirst      Jacobian of the integration operation relative to
         * the state point (size `ndx`\f$\times\f$`ndx`)
         * @param[out] Jsecond     Jacobian of the integration operation relative to
         * the velocity vector (size `ndx`\f$\times\f$`ndx`)
         * @param[in] firstsecond  Argument (either x and / or dx) with respect to
         * which the differentiation is performed
         * @param[in] op           Assignment operator which sets, adds, or removes
         * the given Jacobian matrix
         */
        void Jintegrate(const Eigen::Ref<const VectorXs_t> &x,
                        const Eigen::Ref<const VectorXs_t> &dx,
                        Eigen::Ref<MatrixXs_t> Jfirst,
                        Eigen::Ref<MatrixXs_t> Jsecond,
                        const Jcomponent firstsecond = both,
                        const AssignmentOp op = setto) const
        {
            derived().Jintegrate(x, dx, Jfirst, Jsecond, firstsecond, op);
        }

        /**
         * @brief Parallel transport from integrate(x, dx) to x.
         *
         * This function performs the parallel transportation of an input matrix whose
         * columns are expressed in the tangent space at
         * \f$\mathbf{x}\oplus\delta\mathbf{x}\f$ to the tangent space at
         * \f$\mathbf{x}\f$ point.
         *
         * @param[in]  x           State point (size `nx`).
         * @param[in]  dx          Velocity vector (size `ndx`)
         * @param[out] Jin         Input matrix (number of rows = `nv`).
         * @param[in] firstsecond  Argument (either x or dx) with respect to which the
         * differentiation of Jintegrate is performed.
         */
        void JintegrateTransport(const Eigen::Ref<const VectorXs_t> &x,
                                 const Eigen::Ref<const VectorXs_t> &dx,
                                 Eigen::Ref<MatrixXs_t> Jin,
                                 const Jcomponent firstsecond) const
        {
            derived().JintegrateTransport(x, dx, Jin, firstsecond);
        }

        /**
         * @copybrief diff()
         *
         * @param[in]  x0     Previous state point (size `nx`)
         * @param[in]  x1     Current state point (size `nx`)
         * @return  Difference between the current and previous state points (size
         * `ndx`)
         */
        VectorXs_t diff_dx(const Eigen::Ref<const VectorXs_t> &x0,
                           const Eigen::Ref<const VectorXs_t> &x1)
        {
            derived().diff_dx(x0, x1);
        }

        /**
         * @copybrief integrate()
         *
         * @param[in]  x     State point (size `nx`)
         * @param[in]  dx    Velocity vector (size `ndx`)
         * @return  Next state point (size `nx`)
         */
        VectorXs_t integrate_x(const Eigen::Ref<const VectorXs_t> &x,
                               const Eigen::Ref<const VectorXs_t> &dx)
        {
            derived().integrate_x(x, dx);
        }

        /**
         * @copybrief jdiff()
         *
         * @param[in]  x0     Previous state point (size `nx`)
         * @param[in]  x1     Current state point (size `nx`)
         * @return  Jacobians
         */
        std::vector<MatrixXs_t> Jdiff_Js(const Eigen::Ref<const VectorXs_t> &x0,
                                         const Eigen::Ref<const VectorXs_t> &x1,
                                         const Jcomponent firstsecond = both)
        {
            derived().Jdiff_Js(x0, x1, firstsecond);
        }

        /**
         * @copybrief Jintegrate()
         *
         * @param[in]  x     State point (size `nx`)
         * @param[in]  dx    Velocity vector (size `ndx`)
         * @return  Jacobians
         */
        std::vector<MatrixXs_t> Jintegrate_Js(const Eigen::Ref<const VectorXs_t> &x,
                                              const Eigen::Ref<const VectorXs_t> &dx,
                                              const Jcomponent firstsecond = both)
        {
            derived().Jintegrate_Js(x, dx, firstsecond);
        }

        /**
         * @brief Return the dimension of the state tuple
         */
        std::size_t get_nx() const
        {
            return derived().get_nx();
        }

        /**
         * @brief Return the dimension of the tangent space of the state manifold
         */
        std::size_t get_ndx() const
        {
            return derived().get_ndx();
        }

        /**
         * @brief Return the state lower bound
         */
        const VectorXs_t &get_lb() const
        {
            return derived().get_lb();
        }

        /**
         * @brief Return the state upper bound
         */
        const VectorXs_t &get_ub() const
        {
            return derived().get_ub();
        }

        /**
         * @brief Indicate if the state has defined limits
         */
        bool get_has_limits() const
        {
            return derived().get_has_limits();
        }

        /**
         * @brief Modify the state lower bound
         */
        void set_lb(const VectorXs_t &lb)
        {
            derived().set_lb(lb);
        }

        /**
         * @brief Modify the state upper bound
         */
        void set_ub(const VectorXs_t &ub)
        {
            derived().set_ub(ub);
        }

    protected:
        void update_has_limits()
        {
            derived().update_has_limits();
        }

        // Default constructor: protected.
        // Prevent the construction of stand-alone StateBase.
        inline StateBase() : nx_(0), ndx_(0), has_limits_(false)
        {
            lb_.resize(0);
            ub_.resize(0);
        }

        // Copy constructor: protected.
        // Copy of stand-alone StateBase are prevented, but can be used from inheriting
        // objects. Copy is done by calling copy operator.
        inline StateBase(const StateBase &clone)
        {
            *this = clone;
        }

        // Copy operator: protected.
        // Copy of stand-alone StateBase are prevented, but can be used from inheriting
        // objects.
        inline StateBase &operator=(const StateBase &clone)
        {
            nx_ = clone.nx_;
            ndx_ = clone.ndx_;
            lb_ = clone.lb_;
            ub_ = clone.ub_;
            has_limits_ = clone.has_limits_;
            return *this;
        }

        std::size_t nx_;  //!< State dimension
        std::size_t ndx_; //!< State rate dimension
        VectorXs_t lb_;   //!< Lower state limits
        VectorXs_t ub_;   //!< Upper state limits
        bool has_limits_; //!< Indicates whether any of the state limits is finite
    };

} // namespace galileo