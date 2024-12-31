#pragma once

#include <boost/shared_ptr.hpp>

#include "galileo/fwd.hpp"
#include "galileo/core/math-base.hpp"

#define GALILEO_CONTROL_MODEL_TYPEDEF_GENERIC(Control, TYPENAME)               \
    typedef TYPENAME traits<Control>::Scalar Scalar;                           \
    typedef TYPENAME traits<Control>::ControlModelDerived ControlModelDerived; \
    typedef TYPENAME traits<Control>::ControlDataDerived ControlDataDerived;   \
    typedef TYPENAME traits<Control>::VectorXs_t VectorXs_t;                   \
    typedef TYPENAME traits<Control>::MatrixXs_t MatrixXs_t;

#define GALILEO_CONTROL_TYPEDEF(Control) GALILEO_CONTROL_MODEL_TYPEDEF_GENERIC(Control, typename)
#define GALILEO_CONTROL_TYPEDEF_TEMPLATE(Control) \
    GALILEO_CONTROL_MODEL_TYPEDEF_GENERIC(Control, typename)

#define GALILEO_CONTROL_CAST_TYPE_SPECIALIZATION(ControlModelTpl) \
    template <typename Scalar, typename NewScalar>                \
    struct CastType<NewScalar, ControlModelTpl<Scalar>>           \
    {                                                             \
        typedef ControlModelTpl<NewScalar> type;                  \
    }

namespace galileo
{
    /**
     * @brief Abstract class for the control trajectory parametrization
     *
     * The control trajectory is a function of the (normalized) time.
     * Normalized time is between 0 and 1, where 0 represents the beginning of the
     * time step, and 1 represents its end. The trajectory depends on the control
     * parameters u, whose size may be larger than the size of the control inputs w.
     *
     * The main computations are carried out in `calc`, `multiplyByJacobian` and
     * `multiplyJacobianTransposeBy`, where the former computes control input
     * \f$\mathbf{w}\in\mathbb{R}^{nw}\f$ from a set of control parameters
     * \f$\mathbf{u}\in\mathbb{R}^{nu}\f$ where `nw` and `nu` represent the
     * dimension of the control inputs and parameters, respectively, and the latter
     * defines useful operations across the Jacobian of the control-parametrization
     * model. Finally, `params` allows us to obtain the control parameters from a
     * the control input, i.e., it is the dual of `calc`. Note that
     * `multiplyByJacobian` and `multiplyJacobianTransposeBy` requires to run `calc`
     * first.
     */
    template <typename Derived>
    class ControlModelBase : NumericalBase<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::ControlDerived ControlDerived;
        GALILEO_CONTROL_TYPEDEF_TEMPLATE(ControlDerived);

        ControlModelDerived &derived()
        {
            return *static_cast<Derived *>(this);
        }

        const ControlModelDerived &derived() const
        {
            return *static_cast<const Derived *>(this);
        }

        /**
         * @brief Get the value of the control at the specified time
         *
         * @param[in]  data   Data structure containing the control vector to write
         * @param[in]  t      Time in [0,1]
         * @param[in]  u      Control parameters
         */
        void calc(
            const boost::shared_ptr<ControlParametrizationDataAbstract> &data,
            const Scalar t, const Eigen::Ref<const VectorXs_t> &u) const
        {
            derived().calc(data, t, u);
        }

        /**
         * @brief Get the value of the Jacobian of the control with respect to the
         * parameters
         *
         * It assumes that `calc()` has been run first
         *
         * @param[in]  data   Control-parametrization data
         * @param[in]  t      Time in [0,1]
         * @param[in]  u      Control parameters
         */
        void calcDiff(
            const boost::shared_ptr<ControlParametrizationDataAbstract> &data,
            const Scalar t, const Eigen::Ref<const VectorXs_t> &u) const
        {
            derived().calcDiff(data, t, u);
        }

        /**
         * @brief Create the control-parametrization data
         *
         * @return the control-parametrization data
         */
        boost::shared_ptr<ControlParametrizationDataAbstract> createData()
        {
            return derived().createData();
        }

        /**
         * @brief Update the control parameters u for a specified time t given the
         * control input w
         *
         * @param[in]  data   Control-parametrization data
         * @param[in]  t      Time in [0,1]
         * @param[in]  w      Control inputs
         */
        void params(
            const boost::shared_ptr<ControlParametrizationDataAbstract> &data,
            const Scalar t, const Eigen::Ref<const VectorXs_t> &w) const
        {
            derived().params(data, t, w);
        }

        /**
         * @brief Convert the bounds on the control inputs w to bounds on the control
         * parameters u
         *
         * @param[in]  w_lb   Control lower bound
         * @param[in]  w_ub   Control lower bound
         * @param[out] u_lb   Control parameters lower bound
         * @param[out] u_ub   Control parameters upper bound
         */
        void convertBounds(const Eigen::Ref<const VectorXs_t> &w_lb,
                           const Eigen::Ref<const VectorXs_t> &w_ub,
                           Eigen::Ref<VectorXs_t> u_lb,
                           Eigen::Ref<VectorXs_t> u_ub) const
        {
            derived().convertBounds(w_lb, w_ub, u_lb, u_ub);
        }

        /**
         * @brief Compute the product between the given matrix A and the derivative of
         * the control input with respect to the control parameters (i.e., A*dw_du).
         *
         * It assumes that `calc()` has been run first
         *
         * @param[in]  data   Control-parametrization data
         * @param[in]  A      A matrix to multiply times the Jacobian
         * @param[out] out    Product between the matrix A and the Jacobian of the
         * control with respect to the parameters
         * @param[in] op      Assignment operator which sets, adds, or removes the
         * given results
         */
        void multiplyByJacobian(
            const boost::shared_ptr<ControlParametrizationDataAbstract> &data,
            const Eigen::Ref<const MatrixXs_t> &A, Eigen::Ref<MatrixXs_t> out,
            const AssignmentOp = setto) const
        {
            derived().multiplyByJacobian(data, A, out);
        }

        MatrixXs_t multiplyByJacobian_J(
            const boost::shared_ptr<ControlParametrizationDataAbstract> &data,
            const Eigen::Ref<const MatrixXs_t> &A, const AssignmentOp = setto) const
        {
            return derived().multiplyByJacobian_J(data, A);
        }

        /**
         * @brief Compute the product between the transpose of the derivative of the
         * control input with respect to the control parameters and a given matrix A
         * (i.e., dw_du^T*A)
         *
         * It assumes that `calc()` has been run first
         *
         * @param[in]  data   Control-parametrization data
         * @param[in]  A      A matrix to multiply times the Jacobian
         * @param[out] out    Product between the transposed Jacobian of the control
         * with respect to the parameters and the matrix A
         * @param[in] op      Assignment operator which sets, adds, or removes the
         * given results
         */
        void multiplyJacobianTransposeBy(
            const boost::shared_ptr<ControlParametrizationDataAbstract> &data,
            const Eigen::Ref<const MatrixXs_t> &A, Eigen::Ref<MatrixXs_t> out,
            const AssignmentOp = setto) const
        {
            derived().multiplyJacobianTransposeBy(data, A, out);
        }

        MatrixXs_t multiplyJacobianTransposeBy_J(
            const boost::shared_ptr<ControlParametrizationDataAbstract> &data,
            const Eigen::Ref<const MatrixXs_t> &A, const AssignmentOp = setto) const
        {
            return derived().multiplyJacobianTransposeBy_J(data, A);
        }

        /**
         * @brief Checks that a specific data belongs to this model
         */
        bool checkData(
            const boost::shared_ptr<ControlParametrizationDataAbstract> &data)
        {
            return derived().checkData(data);
        }

        /**
         * @brief Return the dimension of the control inputs
         */
        std::size_t get_nw() const
        {
            return derived().get_nw();
        }

        /**
         * @brief Return the dimension of control parameters
         */
        std::size_t get_nu() const
        {
            return derived().get_nu();
        }

    protected:
        // Default constructor: protected.
        // Prevent the construction of stand-alone ControlModelBase.
        inline ControlModelBase() : nw_(0), nu_(0)
        {
        }

        // Copy constructor: protected.
        // Copy of stand-alone ControlModelBase are prevented, but can be used from inheriting
        // objects. Copy is done by calling copy operator.
        inline ControlModelBase(const ControlModelBase &clone)
        {
            *this = clone;
        }

        // Copy operator: protected.
        // Copy of stand-alone ControlModelBase are prevented, but can be used from inheriting
        // objects.
        inline ControlModelBase &operator=(const ControlModelBase &clone)
        {
            nw_ = clone.nw_;
            nu_ = clone.nu_;
            return *this;
        }

        std::size_t nw_; //!< Control dimension
        std::size_t nu_; //!< Control parameters dimension
    };

    template <typename Derived>
    struct ControlDataBase : NumericalBase<Derived>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::ControlDerived ControlDerived;
        GALILEO_CONTROL_TYPEDEF_TEMPLATE(ControlDerived);

        VectorXs_t w;     //!< value of the differential control
        VectorXs_t u;     //!< value of the control parameters
        MatrixXs_t dw_du; //!< Jacobian of the differential control with respect to the
                          //!< parameters

    protected:
        // Default constructor: protected.
        inline ControlDataBase()
        {
        }
    };

} // namespace galileo