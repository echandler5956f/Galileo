#ifndef __galileo_core_states_state_base_hpp__
#define __galileo_core_states_state_base_hpp__

#include "galileo/core/fwd.hpp"
#include "galileo/multibody/robot-spec.hpp"

namespace galileo
{

    enum Jcomponent
    {
        both = 0,
        first = 1,
        second = 2
    }; // enum Jcomponent

    inline bool is_a_Jcomponent(Jcomponent firstsecond)
    {
        return (firstsecond == first || firstsecond == second || firstsecond == both);
    }

    template <typename Derived, typename RobotSpec>
    class StateBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

        GALILEO_ROBOT_SPEC_EIGEN_TYPES_TYPEDEF(RS);

        /**
         * @brief Generate a zero state
         */
        VectorNx_t zero() const
        {
            return this->derived().zero();
        }

        /**
         * @brief Generate a random state
         */
        VectorNx_t rand() const
        {
            return this->derived().rand();
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
            this->derived().diff(x0.derived(), x1.derived(), dxout.derived());
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
            this->derived().integrate(x.derived(), dx.derived(), xout.derived());
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
            this->derived().Jdiff(x0.derived(), x1.derived(), Jfirst.derived(), Jsecond.derived(), firstsecond);
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
            this->derived().Jintegrate(x.derived(), dx.derived(), Jfirst.derived(), Jsecond.derived(), firstsecond, op);
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
            this->derived().JintegrateTransport(x.derived(), dx.derived(), Jin.derived(), firstsecond);
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
        VectorNdx_t diff_dx(const Eigen::MatrixBase<StateVector1> &x0,
                            const Eigen::MatrixBase<StateVector2> &x1)
        {
            VectorNdx_t dx;
            this->derived().diff(x0.derived(), x1.derived(), dx);
            return dx;
        }

        /**
         * @copybrief integrate()
         *
         * @param[in]  x     State point (size `nx`)
         * @param[in]  dx    Velocity vector (size `ndx`)
         * @return  Next state point (size `nx`)
         */
        template <typename StateVector, typename StateTangentVector>
        VectorNx_t integrate_x(const Eigen::MatrixBase<StateVector> &x,
                               const Eigen::MatrixBase<StateTangentVector> &dx)
        {
            VectorNx_t xout;
            this->derived().integrate(x.derived(), dx.derived(), xout);
            return xout;
        }

        /**
         * @copybrief jdiff()
         *
         * @param[in]  x0     Previous state point (size `nx`)
         * @param[in]  x1     Current state point (size `nx`)
         * @return  Jacobians
         */
        template <typename StateVector1, typename StateVector2>
        std::vector<MatrixNdx_t> Jdiff_Js(const Eigen::MatrixBase<StateVector1> &x0,
                                          const Eigen::MatrixBase<StateVector2> &x1,
                                          const Jcomponent firstsecond = both)
        {
            MatrixNdx_t Jfirst = MatrixNdx_t::Zero(get_ndx(), get_ndx());
            MatrixNdx_t Jsecond = MatrixNdx_t::Zero(get_ndx(), get_ndx());
            std::vector<MatrixNdx_t> Jacs;
            Jdiff(x0, x1, Jfirst, Jsecond, firstsecond);
            switch (firstsecond)
            {
            case both:
                Jacs.push_back(Jfirst);
                Jacs.push_back(Jsecond);
                break;
            case first:
                Jacs.push_back(Jfirst);
                break;
            case second:
                Jacs.push_back(Jsecond);
                break;
            default:
                Jacs.push_back(Jfirst);
                Jacs.push_back(Jsecond);
                break;
            }
            return Jacs;
        }

        /**
         * @copybrief Jintegrate()
         *
         * @param[in]  x     State point (size `nx`)
         * @param[in]  dx    Velocity vector (size `ndx`)
         * @return  Jacobians
         */
        template <typename StateVector, typename StateTangentVector>
        std::vector<MatrixNdx_t> Jintegrate_Js(const Eigen::MatrixBase<StateVector> &x,
                                               const Eigen::MatrixBase<StateTangentVector> &dx,
                                               const Jcomponent firstsecond = both)
        {
            MatrixNdx_t Jfirst = MatrixNdx_t::Zero(get_ndx(), get_ndx());
            MatrixNdx_t Jsecond = MatrixNdx_t::Zero(get_ndx(), get_ndx());
            std::vector<MatrixNdx_t> Jacs;
            Jintegrate(x, dx, Jfirst, Jsecond, firstsecond, setto);
            switch (firstsecond)
            {
            case both:
                Jacs.push_back(Jfirst);
                Jacs.push_back(Jsecond);
                break;
            case first:
                Jacs.push_back(Jfirst);
                break;
            case second:
                Jacs.push_back(Jsecond);
                break;
            default:
                Jacs.push_back(Jfirst);
                Jacs.push_back(Jsecond);
                break;
            }
            return Jacs;
        }

        RS &get_rs()
        {
            return rs_;
        }

        /**
         * @brief Return the dimension (if any) of the floating base configuration space of the state
         */
        constexpr int get_nqb() const
        {
            return this->derived().get_nqb_impl();
        }

        constexpr int get_nqb_impl() const
        {
            if constexpr (RS::DimNQb_t::IsFixed)
            {
                return RS::DimNQb_t::Value;
            }
            else
            {
                return rs_.NQb_dim.value();
            }
        }

        const RS::DimNQb_t &NQbDim() const
        {
            return rs_.NQb_dim;
        }

        /**
         * @brief Return the dimension of the joint configuration space of the state
         */
        constexpr int get_nqj() const
        {
            return this->derived().get_nqj_impl();
        }

        constexpr int get_nqj_impl() const
        {
            if constexpr (RS::DimNQj_t::IsFixed)
            {
                return RS::DimNQj_t::Value;
            }
            else
            {
                return rs_.NQj_dim.value();
            }
        }

        const RS::DimNQj_t &NQjDim() const
        {
            return rs_.NQj_dim;
        }

        /**
         * @brief Return the dimension of the configuration space of the state
         */
        constexpr int get_nq() const
        {
            return this->derived().get_nq_impl();
        }

        constexpr int get_nq_impl() const
        {
            if constexpr (RS::DimNQ_t::IsFixed)
            {
                return RS::DimNQ_t::Value;
            }
            else
            {
                return rs_.NQ_dim.value();
            }
        }

        const RS::DimNQ_t &NQDim() const
        {
            return rs_.NQ_dim;
        }

        /**
         * @brief Return the dimension (if any) of the floating base velocity space of the state
         */
        constexpr int get_nvb() const
        {
            return this->derived().get_nvb_impl();
        }

        constexpr int get_nvb_impl() const
        {
            if constexpr (RS::DimNVb_t::IsFixed)
            {
                return RS::DimNVb_t::Value;
            }
            else
            {
                return rs_.NVb_dim.value();
            }
        }

        const RS::DimNVb_t &NVbDim() const
        {
            return rs_.NVb_dim;
        }

        /**
         * @brief Return the dimension of the joint velocity space of the state
         */
        constexpr int get_nvj() const
        {
            return this->derived().get_nvj_impl();
        }

        constexpr int get_nvj_impl() const
        {
            if constexpr (RS::DimNVj_t::IsFixed)
            {
                return RS::DimNVj_t::Value;
            }
            else
            {
                return rs_.NVj_dim.value();
            }
        }

        const RS::DimNVj_t &NVjDim() const
        {
            return rs_.NVj_dim;
        }

        /**
         * @brief Return the dimension of the velocity space of the state
         */
        constexpr int get_nv() const
        {
            return this->derived().get_nv_impl();
        }

        constexpr int get_nv_impl() const
        {
            if constexpr (RS::DimNV_t::IsFixed)
            {
                return RS::DimNV_t::Value;
            }
            else
            {
                return rs_.NV_dim.value();
            }
        }

        const RS::DimNV_t &NVDim() const
        {
            return rs_.NV_dim;
        }

        /**
         * @brief Return the number of rotors attached to the floating base (if any)
         */
        constexpr int get_nrotors() const
        {
            return this->derived().get_nrotors_impl();
        }

        constexpr int get_nrotors_impl() const
        {
            if constexpr (RS::DimNRotors_t::IsFixed)
            {
                return RS::DimNRotors_t::Value;
            }
            else
            {
                return rs_.NRotors_dim.value();
            }
        }

        const RS::DimNRotors_t &NRotorsDim() const
        {
            return rs_.NRotors_dim;
        }

        /**
         * @brief Return the dimension of the state
         */
        constexpr int get_nx() const
        {
            return this->derived().get_nx_impl();
        }

        constexpr int get_nx_impl() const
        {
            if constexpr (RS::DimNX_t::IsFixed)
            {
                return RS::DimNX_t::Value;
            }
            else
            {
                return rs_.NX_dim.value();
            }
        }

        const RS::DimNX_t &NXDim() const
        {
            return rs_.NX_dim;
        }

        /**
         * @brief Return the dimension of the tangent space of the state manifold
         */
        constexpr int get_ndx() const
        {
            return this->derived().get_ndx_impl();
        }

        constexpr int get_ndx_impl() const
        {
            if constexpr (RS::DimNDX_t::IsFixed)
            {
                return RS::DimNDX_t::Value;
            }
            else
            {
                return rs_.NDX_dim.value();
            }
        }

        const RS::DimNDX_t &NDXDim() const
        {
            return rs_.NDX_dim;
        }

        /**
         * @brief Return the dimension of the actuated torque space of the state
         */
        constexpr int get_nua() const
        {
            return this->derived().get_nua_impl();
        }

        constexpr int get_nua_impl() const
        {
            if constexpr (RS::DimNUa_t::IsFixed)
            {
                return RS::DimNUa_t::Value;
            }
            else
            {
                return rs_.NUa_dim.value();
            }
        }

        const RS::DimNUa_t &NUaDim() const
        {
            return rs_.NUa_dim;
        }

        /**
         * @brief Return the state lower bound
         */
        const VectorNx_t &get_lb() const
        {
            return this->derived().get_lb_impl();
        }

        /**
         * @brief Return the state upper bound
         */
        const VectorNx_t &get_ub() const
        {
            return this->derived().get_ub_impl();
        }

        /**
         * @brief Modify the state lower bound
         */
        template <typename StateVector>
        void set_lb(const Eigen::MatrixBase<StateVector> &lb)
        {
            this->derived().set_lb_impl(lb.derived());
        }

        /**
         * @brief Modify the state upper bound
         */
        template <typename StateVector>
        void set_ub(const Eigen::MatrixBase<StateVector> &ub)
        {
            this->derived().set_ub_impl(ub.derived());
        }

    protected:
        inline StateBase(const RS &rs)
            : rs_(rs)
        {
        }

        inline StateBase(const StateBase &clone)
        {
            *this = clone;
        }

        inline StateBase &operator=(const StateBase &clone)
        {
            return *this;
        }

        RS rs_;

    }; // class StateBase

} // namespace galileo

#endif // __galileo_core_states_state_base_hpp__
