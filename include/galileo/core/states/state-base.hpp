#ifndef __galileo_core_states_state_base_hpp__
#define __galileo_core_states_state_base_hpp__

#include "galileo/core/fwd.hpp"
#include "galileo/multibody/robot-spec.hpp"

namespace galileo
{

    template <typename Derived, typename RobotSpec>
    class StateBase : public internal::CRTP<Derived>
    {
    public:
        using RS = RobotSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(RS);

        /**
         * @brief Generate a zero state
         */
        VectorNx_t zero() const { return this->derived().zero(); }

        /**
         * @brief Generate a random state
         */
        VectorNx_t rand() const { return this->derived().rand(); }

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
            this->derived().diff(x0, x1, dxout);
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
            this->derived().integrate(x, dx, xout);
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
         */
        template <Jcomponent jc, typename StateVector1, typename StateVector2, typename JMatrix1, typename JMatrix2>
        void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                   const Eigen::MatrixBase<StateVector2> &x1,
                   Eigen::MatrixBase<JMatrix1> &Jfirst,
                   Eigen::MatrixBase<JMatrix2> &Jsecond) const
        {
            this->derived().template Jdiff<jc>(x0, x1, Jfirst, Jsecond);
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
         */
        template <Jcomponent jc = BOTH,
                  AssignmentOp op = SETTO,
                  typename StateVector,
                  typename StateTangentVector,
                  typename JMatrix1,
                  typename JMatrix2>
        void Jintegrate(const Eigen::MatrixBase<StateVector> &x,
                        const Eigen::MatrixBase<StateTangentVector> &dx,
                        Eigen::MatrixBase<JMatrix1> &Jfirst,
                        Eigen::MatrixBase<JMatrix2> &Jsecond) const
        {
            this->derived().template Jintegrate<jc, op>(x, dx, Jfirst, Jsecond);
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
         */
        template <Jcomponent jc, typename StateVector, typename StateTangentVector, typename JMatrix>
        void JintegrateTransport(const Eigen::MatrixBase<StateVector> &x,
                                 const Eigen::MatrixBase<StateTangentVector> &dx,
                                 Eigen::MatrixBase<JMatrix> &Jin) const
        {
            this->derived().template JintegrateTransport<jc>(x, dx, Jin);
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
        VectorNdx_t diff_dx(const Eigen::MatrixBase<StateVector1> &x0, const Eigen::MatrixBase<StateVector2> &x1) const
        {
            VectorNdx_t dx = VectorNdx_t::Zero(get_ndx());
            this->derived().diff(x0, x1, dx);
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
                               const Eigen::MatrixBase<StateTangentVector> &dx) const
        {
            VectorNx_t xout = VectorNx_t::Zero(get_nx());
            this->derived().integrate(x, dx, xout);
            return xout;
        }

        /**
         * @copybrief jdiff()
         *
         * @param[in]  x0     Previous state point (size `nx`)
         * @param[in]  x1     Current state point (size `nx`)
         * @return  Jacobians
         */
        template <Jcomponent jc = BOTH, typename StateVector1, typename StateVector2>
        std::vector<MatrixNdx_t> Jdiff_Js(const Eigen::MatrixBase<StateVector1> &x0,
                                          const Eigen::MatrixBase<StateVector2> &x1) const
        {
            MatrixNdx_t Jfirst = MatrixNdx_t::Zero(get_ndx(), get_ndx());
            MatrixNdx_t Jsecond = MatrixNdx_t::Zero(get_ndx(), get_ndx());
            std::vector<MatrixNdx_t> Jacs;
            Jdiff<jc>(x0, x1, Jfirst, Jsecond);

            if constexpr (IsFirst<jc> || IsBoth<jc>) Jacs.push_back(Jfirst);
            if constexpr (IsSecond<jc> || IsBoth<jc>) Jacs.push_back(Jsecond);
            return Jacs;
        }

        /**
         * @copybrief Jintegrate()
         *
         * @param[in]  x     State point (size `nx`)
         * @param[in]  dx    Velocity vector (size `ndx`)
         * @return  Jacobians
         */
        template <Jcomponent jc = BOTH, typename StateVector, typename StateTangentVector>
        std::vector<MatrixNdx_t> Jintegrate_Js(const Eigen::MatrixBase<StateVector> &x,
                                               const Eigen::MatrixBase<StateTangentVector> &dx) const
        {
            MatrixNdx_t Jfirst = MatrixNdx_t::Zero(get_ndx(), get_ndx());
            MatrixNdx_t Jsecond = MatrixNdx_t::Zero(get_ndx(), get_ndx());
            std::vector<MatrixNdx_t> Jacs;
            Jintegrate<jc, SETTO>(x, dx, Jfirst, Jsecond);

            if constexpr (IsFirst<jc> || IsBoth<jc>) Jacs.push_back(Jfirst);
            if constexpr (IsSecond<jc> || IsBoth<jc>) Jacs.push_back(Jsecond);

            return Jacs;
        }

        const RS &get_rs() const { return rs_; }
        RS &get_rs() { return rs_; }

        /**
         * @brief Return the dimension (if any) of the floating base configuration space of the state
         */
        const DimNQb_t &get_nqb_dim() const { return rs_.get_nqb_dim(); }
        int get_nqb() const { return rs_.get_nqb(); }

        /**
         * @brief Return the dimension of the joint configuration space of the state
         */
        const DimNQj_t &get_nqj_dim() const { return rs_.get_nqj_dim(); }
        int get_nqj() const { return rs_.get_nqj(); }

        /**
         * @brief Return the dimension of the configuration space of the state
         */
        const DimNQ_t &get_nq_dim() const { return rs_.get_nq_dim(); }
        int get_nq() const { return rs_.get_nq(); }

        /**
         * @brief Return the dimension (if any) of the floating base velocity space of the state
         */
        const DimNVb_t &get_nvb_dim() const { return rs_.get_nvb_dim(); }
        int get_nvb() const { return rs_.get_nvb(); }

        /**
         * @brief Return the dimension of the joint velocity space of the state
         */
        const DimNVj_t &get_nvj_dim() const { return rs_.get_nvj_dim(); }
        int get_nvj() const { return rs_.get_nvj(); }

        /**
         * @brief Return the dimension of the velocity space of the state
         */
        const DimNV_t &get_nv_dim() const { return rs_.get_nv_dim(); }
        int get_nv() const { return rs_.get_nv(); }

        /**
         * @brief Return the number of rotors attached to the floating base (if any)
         */
        const DimNRotors_t &get_nrotors_dim() const { return rs_.get_nrotors_dim(); }
        int get_nrotors() const { return rs_.get_nrotors(); }

        /**
         * @brief Return the dimension of the state
         */
        const DimNX_t &get_nx_dim() const { return rs_.get_nx_dim(); }
        int get_nx() const { return rs_.get_nx(); }

        /**
         * @brief Return the dimension of the tangent space of the state manifold
         */
        const DimNDX_t &get_ndx_dim() const { return rs_.get_ndx_dim(); }
        int get_ndx() const { return rs_.get_ndx(); }

        /**
         * @brief Return the dimension of the actuated torque space of the state
         */
        const DimNUa_t &get_nua_dim() const { return rs_.get_nua_dim(); }
        int get_nua() const { return rs_.get_nua(); }

        /**
         * @brief Return the state lower bound
         */
        const VectorNx_t &get_lb() const { return this->derived().get_lb(); }

        /**
         * @brief Return the state upper bound
         */
        const VectorNx_t &get_ub() const { return this->derived().get_ub(); }

        /**
         * @brief Modify the state lower bound
         */
        template <typename StateVector>
        void set_lb(const Eigen::MatrixBase<StateVector> &lb)
        {
            this->derived().set_lb(lb);
        }

        /**
         * @brief Modify the state upper bound
         */
        template <typename StateVector>
        void set_ub(const Eigen::MatrixBase<StateVector> &ub)
        {
            this->derived().set_ub(ub);
        }

        /**
         * @brief Display state information to output stream
         * @param os Output stream
         *
         * This method must be implemented by derived classes to provide
         * state-specific display information.
         */
        void display(std::ostream &os, const std::string &indent = "  ") const { this->derived().display(os, indent); }

        // Stream output operator
        friend std::ostream &operator<<(std::ostream &os, const StateBase &state)
        {
            os << "State: {\n";
            state.derived().display(os, "  ");
            os << "}";
            return os;
        }

    protected:
        inline StateBase(const RS &rs) : rs_(rs) {}
        inline StateBase(const StateBase &clone) : rs_(clone.rs_) {}
        inline StateBase &operator=(const StateBase &clone)
        {
            rs_ = clone.rs_;
            return *this;
        }

        RS rs_;

    }; // class StateBase

} // namespace galileo

#endif // __galileo_core_states_state_base_hpp__
