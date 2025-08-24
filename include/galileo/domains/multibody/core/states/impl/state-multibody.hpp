#ifndef __galileo_multibody_core_states_state_multibody_hpp__
#define __galileo_multibody_core_states_state_multibody_hpp__

#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/multibody/model.hpp>

#include "galileo/core/states/state-base.hpp"
#include "galileo/domains/multibody/core/states/fwd.hpp"

#include <limits>

namespace galileo
{

    template <typename MultibodySpec>
    class StateMultibodyTpl : public StateBase<StateMultibodyTpl<MultibodySpec>, MultibodySpec>
    {
    public:
        using SS = MultibodySpec;

        using Base = StateBase<StateMultibodyTpl<SS>, SS>;

        using RobotModel_t = typename SS::RobotModel_t;
        using VectorNx_t = typename SS::VectorNx_t;
        using VectorNqb_t = typename SS::VectorNqb_t;
        using VarScalar = typename SS::VarScalar;
        using NumScalar = typename SS::VarScalar;

        StateMultibodyTpl(const SS &ss, const RobotModel_t &model) : Base(ss), model_(model) { initialize(); }

        VectorNx_t zero() const { return x0_; }

        VectorNx_t rand() const
        {
            VectorNx_t xrand = VectorNx_t::Random(get_ss().get_nx());
            head(xrand, get_ss().get_nq_dim()) = pinocchio::randomConfiguration(model_);

            // For the 3x1 position component, set to a uniform random distribution
            // between -1 and 1
            // Need to add checks based on the type of the first joint.
            // Currently assumes pinocchio::JointModelFreeFlyer
            // TODO: Add support for other joint types
            if (get_ss().get_nqb() >= 3)
            {
                head(xrand, 3) = Eigen::Matrix<NumScalar, 3, 1>::Random();
            }

            return xrand;
        }

        template <typename StateVector1, typename StateVector2, typename StateTangentVector>
        void diff(const Eigen::MatrixBase<StateVector1> &x0,
                  const Eigen::MatrixBase<StateVector2> &x1,
                  Eigen::MatrixBase<StateTangentVector> &dxout) const
        {
            pinocchio::difference(model_,
                                  head(x0, get_ss().get_nq_dim()),
                                  head(x1, get_ss().get_nq_dim()),
                                  head(dxout, get_ss().get_nv_dim()));
            tail(dxout, get_ss().get_nv_dim()) = tail(x1, get_ss().get_nv_dim()) - tail(x0, get_ss().get_nv_dim());
        }

        template <typename StateVector1, typename StateTangentVector, typename StateVector2>
        void integrate(const Eigen::MatrixBase<StateVector1> &x,
                       const Eigen::MatrixBase<StateTangentVector> &dx,
                       Eigen::MatrixBase<StateVector2> &xout) const
        {
            pinocchio::integrate(model_,
                                 head(x, get_ss().get_nq_dim()),
                                 head(dx, get_ss().get_nv_dim()),
                                 head(xout, get_ss().get_nq_dim()));
            tail(xout, get_ss().get_nv_dim()) = tail(x, get_ss().get_nv_dim()) + tail(dx, get_ss().get_nv_dim());
        }

        template <Jcomponent jc = BOTH,
                  typename StateVector1,
                  typename StateVector2,
                  typename JMatrix1,
                  typename JMatrix2>
        void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                   const Eigen::MatrixBase<StateVector2> &x1,
                   Eigen::MatrixBase<JMatrix1> &Jfirst,
                   Eigen::MatrixBase<JMatrix2> &Jsecond) const
        {
            if constexpr (IsFirst<jc> || IsBoth<jc>) Jfirst.setZero();
            if constexpr (IsSecond<jc> || IsBoth<jc>) Jsecond.setZero();

            if constexpr (IsFirst<jc> || IsBoth<jc>)
            {
                pinocchio::dDifference(model_,
                                       head(x0, get_ss().get_nq_dim()),
                                       head(x1, get_ss().get_nq_dim()),
                                       topLeftCorner(Jfirst, get_ss().get_nv_dim(), get_ss().get_nv_dim()),
                                       pinocchio::ARG0);
                bottomRightCorner(Jfirst, get_ss().get_nv_dim(), get_ss().get_nv_dim()).diagonal().array() =
                    VarScalar(-1.);
            }
            if constexpr (IsSecond<jc> || IsBoth<jc>)
            {
                pinocchio::dDifference(model_,
                                       head(x0, get_ss().get_nq_dim()),
                                       head(x1, get_ss().get_nq_dim()),
                                       topLeftCorner(Jsecond, get_ss().get_nv_dim(), get_ss().get_nv_dim()),
                                       pinocchio::ARG1);
                bottomRightCorner(Jsecond, get_ss().get_nv_dim(), get_ss().get_nv_dim()).diagonal().array() =
                    VarScalar(1.);
            }
        }

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
            // Only zero the matrices for SETTO operations, not for ADDTO/RMFROM
            if constexpr (IsSetTo<op>)
            {
                if constexpr (IsFirst<jc> || IsBoth<jc>) Jfirst.setZero();
                if constexpr (IsSecond<jc> || IsBoth<jc>) Jsecond.setZero();
            }

            if constexpr (IsFirst<jc> || IsBoth<jc>)
            {
                if constexpr (IsSetTo<op>)
                {
                    pinocchio::dIntegrate(model_,
                                          head(x, get_ss().get_nq_dim()),
                                          head(dx, get_ss().get_nv_dim()),
                                          topLeftCorner(Jfirst, get_ss().get_nv_dim(), get_ss().get_nv_dim()),
                                          pinocchio::ARG0,
                                          pinocchio::SETTO);
                    bottomRightCorner(Jfirst, get_ss().get_nv_dim(), get_ss().get_nv_dim()).diagonal().array() =
                        VarScalar(1.);
                }
                else if constexpr (IsAddTo<op>)
                {
                    pinocchio::dIntegrate(model_,
                                          head(x, get_ss().get_nq_dim()),
                                          head(dx, get_ss().get_nv_dim()),
                                          topLeftCorner(Jfirst, get_ss().get_nv_dim(), get_ss().get_nv_dim()),
                                          pinocchio::ARG0,
                                          pinocchio::ADDTO);
                    bottomRightCorner(Jfirst, get_ss().get_nv_dim(), get_ss().get_nv_dim()).diagonal().array() +=
                        VarScalar(1.);
                }
                else if constexpr (IsRmFrom<op>)
                {
                    pinocchio::dIntegrate(model_,
                                          head(x, get_ss().get_nq_dim()),
                                          head(dx, get_ss().get_nv_dim()),
                                          topLeftCorner(Jfirst, get_ss().get_nv_dim(), get_ss().get_nv_dim()),
                                          pinocchio::ARG0,
                                          pinocchio::RMTO);
                    bottomRightCorner(Jfirst, get_ss().get_nv_dim(), get_ss().get_nv_dim()).diagonal().array() -=
                        VarScalar(1.);
                }
            }
            if constexpr (IsSecond<jc> || IsBoth<jc>)
            {
                if constexpr (IsSetTo<op>)
                {
                    pinocchio::dIntegrate(model_,
                                          head(x, get_ss().get_nq_dim()),
                                          head(dx, get_ss().get_nv_dim()),
                                          topLeftCorner(Jsecond, get_ss().get_nv_dim(), get_ss().get_nv_dim()),
                                          pinocchio::ARG1,
                                          pinocchio::SETTO);
                    bottomRightCorner(Jsecond, get_ss().get_nv_dim(), get_ss().get_nv_dim()).diagonal().array() =
                        VarScalar(1.);
                }
                else if constexpr (IsAddTo<op>)
                {
                    pinocchio::dIntegrate(model_,
                                          head(x, get_ss().get_nq_dim()),
                                          head(dx, get_ss().get_nv_dim()),
                                          topLeftCorner(Jsecond, get_ss().get_nv_dim(), get_ss().get_nv_dim()),
                                          pinocchio::ARG1,
                                          pinocchio::ADDTO);
                    bottomRightCorner(Jsecond, get_ss().get_nv_dim(), get_ss().get_nv_dim()).diagonal().array() +=
                        VarScalar(1.);
                }
                else if constexpr (IsRmFrom<op>)
                {
                    pinocchio::dIntegrate(model_,
                                          head(x, get_ss().get_nq_dim()),
                                          head(dx, get_ss().get_nv_dim()),
                                          topLeftCorner(Jsecond, get_ss().get_nv_dim(), get_ss().get_nv_dim()),
                                          pinocchio::ARG1,
                                          pinocchio::RMTO);
                    bottomRightCorner(Jsecond, get_ss().get_nv_dim(), get_ss().get_nv_dim()).diagonal().array() -=
                        VarScalar(1.);
                }
            }
        }

        template <Jcomponent jc, typename StateVector, typename StateTangentVector, typename JMatrix>
        void JintegrateTransport(const Eigen::MatrixBase<StateVector> &x,
                                 const Eigen::MatrixBase<StateTangentVector> &dx,
                                 Eigen::MatrixBase<JMatrix> &Jin) const
        {
            if constexpr (IsFirst<jc> || IsBoth<jc>)
            {
                pinocchio::dIntegrateTransport(model_,
                                               head(x, get_ss().get_nq_dim()),
                                               head(dx, get_ss().get_nv_dim()),
                                               topRows(Jin, get_ss().get_nv_dim()),
                                               pinocchio::ARG0);
            }
            if constexpr (IsSecond<jc> || IsBoth<jc>)
            {
                pinocchio::dIntegrateTransport(model_,
                                               head(x, get_ss().get_nq_dim()),
                                               head(dx, get_ss().get_nv_dim()),
                                               topRows(Jin, get_ss().get_nv_dim()),
                                               pinocchio::ARG1);
            }
        }

        const RobotModel_t &get_robot() const { return model_; }
        const VectorNx_t &get_x0() const { return x0_; }
        const VectorNx_t &get_lb() const { return lb_; }
        const VectorNx_t &get_ub() const { return ub_; }

        template <typename StateVector>
        void set_lb(const Eigen::MatrixBase<StateVector> &lb)
        {
            lb_ = lb;
        }

        template <typename StateVector>
        void set_ub(const Eigen::MatrixBase<StateVector> &ub)
        {
            ub_ = ub;
        }

        using Base::diff_dx;
        using Base::integrate_x;
        using Base::Jdiff_Js;
        using Base::Jintegrate_Js;
        using Base::get_ss;

        /**
         * @brief Display multibody state information to output stream
         */
        void display(std::ostream &os, const std::string &indent = "  ") const
        {
            os << indent << "MultibodySpec: {\n";
            get_ss().display(os, indent + "  ");
            os << indent << "}\n";
            os << indent << "Pinocchio Model: {\n";
            os << indent << "  Model name: " << model_.name << "\n";
            os << indent << "  Number of joints: " << model_.njoints << "\n";
            os << indent << "  Number of bodies: " << model_.nbodies << "\n";
            os << indent << "  Number of frames: " << model_.nframes << "\n";
            os << indent << "  Has floating base: " << (get_ss().get_nqb() > 0 ? "Yes" : "No") << "\n";
            os << indent << "}\n";
            os << indent << "State Bounds: {\n";
            os << indent << "  Lower bounds: " << lb_.transpose() << "\n";
            os << indent << "  Upper bounds: " << ub_.transpose() << "\n";
            os << indent << "  Default state: " << x0_.transpose() << "\n";
            os << indent << "}\n";
        }

    protected:
        void initialize()
        {
            // Now that all of the dynamic dimensions have been updated in the constructor,
            // we can initialize the member variables
            x0_ = VectorNx_t::Zero(get_ss().get_nx());
            head(x0_, get_ss().get_nq_dim()) = pinocchio::neutral(model_);
            lb_ = -VectorNx_t::Constant(get_ss().get_nx(), std::numeric_limits<NumScalar>::infinity());
            ub_ = VectorNx_t::Constant(get_ss().get_nx(), std::numeric_limits<NumScalar>::infinity());

            head(lb_, get_ss().get_nqb_dim()) =
                -VectorNqb_t::Constant(get_ss().get_nqb(), std::numeric_limits<NumScalar>::max());
            head(ub_, get_ss().get_nqb_dim()) =
                VectorNqb_t::Constant(get_ss().get_nqb(), std::numeric_limits<NumScalar>::max());

            segment(lb_, get_ss().get_nqb(), get_ss().get_nqj_dim()) =
                tail(model_.lowerPositionLimit, get_ss().get_nqj_dim());
            segment(ub_, get_ss().get_nqb(), get_ss().get_nqj_dim()) =
                tail(model_.upperPositionLimit, get_ss().get_nqj_dim());

            segment(lb_, get_ss().get_nq(), get_ss().get_nv_dim()) = -model_.velocityLimit;
            segment(ub_, get_ss().get_nq(), get_ss().get_nv_dim()) = model_.velocityLimit;
        }

        RobotModel_t model_;
        VectorNx_t x0_;
        VectorNx_t lb_;
        VectorNx_t ub_;

    }; // class StateMultibodyTpl

} // namespace galileo

#endif // __galileo_multibody_core_states_state_multibody_hpp__
