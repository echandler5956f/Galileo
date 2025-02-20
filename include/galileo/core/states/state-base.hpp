#ifndef __galileo_core_states_state_base_hpp__
#define __galileo_core_states_state_base_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    namespace core
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

        enum AssignmentOp
        {
            setto,
            addto,
            rmfrom
        };

        inline bool is_a_AssignmentOp(AssignmentOp op)
        {
            return (op == setto || op == addto || op == rmfrom);
        }

        template <typename Derived>
        class StateBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
            template <typename StateVector1, typename StateVector2, typename StateTangentVector>
            void diff(const Eigen::MatrixBase<StateVector1> &x0,
                      const Eigen::MatrixBase<StateVector2> &x1,
                      Eigen::MatrixBase<StateTangentVector> &dxout) const
            {
                derived().diff(x0.derived(), x1.derived(), dxout.derived());
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
            template <typename StateVector1, typename StateTangentVector, typename StateVector2>
            void integrate(const Eigen::MatrixBase<StateVector1> &x,
                           const Eigen::MatrixBase<StateTangentVector> &dx,
                           Eigen::MatrixBase<StateVector2> &xout) const
            {
                derived().integrate(x.derived(), dx.derived(), xout.derived());
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
            template <typename StateVector1, typename StateVector2, typename JMatrix1, typename JMatrix2>
            void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                       const Eigen::MatrixBase<StateVector2> &x1,
                       Eigen::MatrixBase<JMatrix1> &Jfirst, Eigen::MatrixBase<JMatrix2> &Jsecond,
                       const Jcomponent firstsecond = both) const
            {
                derived().Jdiff(x0.derived(), x1.derived(), Jfirst.derived(), Jsecond.derived(), firstsecond);
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
            template <typename StateVector, typename StateTangentVector, typename JMatrix1, typename JMatrix2>
            void Jintegrate(const Eigen::MatrixBase<StateVector> &x,
                            const Eigen::MatrixBase<StateTangentVector> &dx,
                            Eigen::MatrixBase<JMatrix1> &Jfirst,
                            Eigen::MatrixBase<JMatrix2> &Jsecond,
                            const Jcomponent firstsecond = both,
                            const AssignmentOp op = setto) const
            {
                derived().Jintegrate(x.derived(), dx.derived(), Jfirst.derived(), Jsecond.derived(), firstsecond, op);
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
             * @param[out] Jin         Input matrix (number of rows = `nv`, number of columns = `ndx`)
             * @param[in] firstsecond  Argument (either x or dx) with respect to which the
             * differentiation of Jintegrate is performed.
             */
            template <typename StateVector, typename StateTangentVector, typename JMatrix>
            void JintegrateTransport(const Eigen::MatrixBase<StateVector> &x,
                                     const Eigen::MatrixBase<StateTangentVector> &dx,
                                     Eigen::MatrixBase<JMatrix> &Jin,
                                     const Jcomponent firstsecond) const
            {
                derived().JintegrateTransport(x.derived(), dx.derived(), Jin.derived(), firstsecond);
            }

            /**
             * @copybrief diff()
             *
             * @param[in]  x0     Previous state point (size `nx`)
             * @param[in]  x1     Current state point (size `nx`)
             * @return  Difference between the current and previous state points (size
             * `ndx`)
             */
            template <typename StateVector1, typename StateVector2>
            VectorXs_t diff_dx(const Eigen::MatrixBase<StateVector1> &x0,
                               const Eigen::MatrixBase<StateVector2> &x1)
            {
                derived().diff_dx(x0.derived(), x1.derived());
            }

            /**
             * @copybrief integrate()
             *
             * @param[in]  x     State point (size `nx`)
             * @param[in]  dx    Velocity vector (size `ndx`)
             * @return  Next state point (size `nx`)
             */
            template <typename StateVector, typename StateTangentVector>
            VectorXs_t integrate_x(const Eigen::MatrixBase<StateVector> &x,
                                   const Eigen::MatrixBase<StateTangentVector> &dx)
            {
                derived().integrate_x(x.derived(), dx.derived());
            }

            /**
             * @copybrief jdiff()
             *
             * @param[in]  x0     Previous state point (size `nx`)
             * @param[in]  x1     Current state point (size `nx`)
             * @return  Jacobians
             */
            template <typename StateVector1, typename StateVector2>
            std::vector<MatrixXs_t> Jdiff_Js(const Eigen::MatrixBase<StateVector1> &x0,
                                             const Eigen::MatrixBase<StateVector2> &x1,
                                             const Jcomponent firstsecond = both)
            {
                derived().Jdiff_Js(x0.derived(), x1.derived(), firstsecond);
            }

            /**
             * @copybrief Jintegrate()
             *
             * @param[in]  x     State point (size `nx`)
             * @param[in]  dx    Velocity vector (size `ndx`)
             * @return  Jacobians
             */
            template <typename StateVector, typename StateTangentVector>
            std::vector<MatrixXs_t> Jintegrate_Js(const Eigen::MatrixBase<StateVector> &x,
                                                  const Eigen::MatrixBase<StateTangentVector> &dx,
                                                  const Jcomponent firstsecond = both)
            {
                derived().Jintegrate_Js(x.derived(), dx.derived(), firstsecond);
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
             * @brief Modify the state lower bound
             */
            template <typename StateVector>
            void set_lb(const Eigen::MatrixBase<StateVector> &lb)
            {
                derived().set_lb(lb.derived());
            }

            /**
             * @brief Modify the state upper bound
             */
            template <typename StateVector>
            void set_ub(const Eigen::MatrixBase<StateVector> &ub)
            {
                derived().set_ub(ub.derived());
            }

        protected:
            inline StateBase() : nx_(0), ndx_(0), has_limits_(false)
            {
            }

            inline StateBase(const StateBase &clone)
            {
                *this = clone;
            }

            inline StateBase &operator=(const StateBase &clone)
            {
                nx_ = clone.nx_;
                ndx_ = clone.ndx_;
                lb_ = clone.lb_;
                ub_ = clone.ub_;
                return *this;
            }

            std::size_t nx_;  //!< State dimension
            std::size_t ndx_; //!< State rate dimension
            VectorXs_t lb_;   //!< Lower state limits
            VectorXs_t ub_;   //!< Upper state limits

        }; // class StateBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_states_state_base_hpp__