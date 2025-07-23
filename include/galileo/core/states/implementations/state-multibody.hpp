#ifndef __galileo_core_states_state_multibody_hpp__
#define __galileo_core_states_state_multibody_hpp__

#include "galileo/core/states/state-base.hpp"
#include "galileo/multibody/robot-spec.hpp"

#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/multibody/model.hpp>

#include <cmath>  // For std::isfinite
#include <limits> // For std::numeric_limits

namespace galileo
{

    template <typename RobotSpec>
    class StateMultibodyTpl
        : public StateBase<StateMultibodyTpl<RobotSpec>, RobotSpec>
    {
    public:
        using RS = RobotSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(RS);

        using Base = StateBase<StateMultibodyTpl<RS>, RS>;

        StateMultibodyTpl(const RobotModel_t &model)
            : Base(RS()),
              model_(model)
        {
            if constexpr (RS::DimNQ_t::IsDynamic)
            {
                get_rs().nq_dim_.set_value(model_.nq);
            }
            if constexpr (RS::DimNQb_t::IsDynamic)
            {
                get_rs().nqb_dim_.set_value(model_.existJointName("root_joint")
                                                ? model_.joints[model_.getJointId("root_joint")].nq()
                                                : 0);
            }
            if constexpr (RS::DimNQj_t::IsDynamic)
            {
                get_rs().nqj_dim_.set_value(get_nq() - get_nqb());
            }
            if constexpr (RS::DimNV_t::IsDynamic)
            {
                get_rs().nv_dim_.set_value(model_.nv);
            }
            if constexpr (RS::DimNVb_t::IsDynamic)
            {
                const std::size_t nvb =
                    model_.existJointName("root_joint")
                        ? model_.joints[model_.getJointId("root_joint")].nv()
                        : 0;
                get_rs().nvb_dim_.set_value(nvb);
            }
            if constexpr (RS::DimNVj_t::IsDynamic)
            {
                const std::size_t nvj = get_nv() - get_nvb();
                get_rs().nvj_dim_.set_value(nvj);
            }
            if constexpr (RS::DimNX_t::IsDynamic)
            {
                get_rs().nx_dim_.set_value(get_nq_dim() + get_nv_dim());
            }
            if constexpr (RS::DimNDX_t::IsDynamic)
            {
                get_rs().ndx_dim_.set_value(get_nv_dim() + get_nv_dim());
            }
            if constexpr (RS::DimNRotors_t::IsDynamic)
            {
                // We do not know anything about rotors at this point in the OCP data pipeline, so if rotors are
                // used, any dynamic evaluation is deferred to ActuationModelFloatingBaseThrustersTpl's constructor.
                // Of course, compile-time NRotors will always be available at this point.
                get_rs().nrotors_dim_.set_value(0);
            }
            if constexpr (RS::DimNUa_t::IsDynamic)
            {
                get_rs().nua_dim_.set_value(get_nvj_dim() + get_nrotors_dim());
            }

            GALILEO_ASSERT(IsValidRobotSpec(get_rs()), "StateMultibodyTpl: Invalid robot spec");
            initialize();
        }

        VectorNx_t zero() const
        {
            return x0_;
        }

        VectorNx_t rand() const
        {
            VectorNx_t xrand = VectorNx_t::Random(get_nx());
            head(xrand, get_nq_dim()) = pinocchio::randomConfiguration(model_);

            // For the 3x1 position component, set to a uniform random distribution
            // between -1 and 1
            // Need to add checks based on the type of the first joint.
            // Currently assumes pinocchio::JointModelFreeFlyer
            // TODO: Add support for other joint types
            if (get_nqb() >= 3)
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
            pinocchio::difference(model_, head(x0, get_nq_dim()), head(x1, get_nq_dim()),
                                  head(dxout, get_nv_dim()));
            tail(dxout, get_nv_dim()) = tail(x1, get_nv_dim()) - tail(x0, get_nv_dim());
        }

        template <typename StateVector1, typename StateTangentVector, typename StateVector2>
        void integrate(const Eigen::MatrixBase<StateVector1> &x,
                       const Eigen::MatrixBase<StateTangentVector> &dx,
                       Eigen::MatrixBase<StateVector2> &xout) const
        {
            pinocchio::integrate(model_, head(x, get_nq_dim()), head(dx, get_nv_dim()),
                                 head(xout, get_nq_dim()));
            tail(xout, get_nv_dim()) = tail(x, get_nv_dim()) + tail(dx, get_nv_dim());
        }

        template <typename StateVector1, typename StateVector2, typename JMatrix1, typename JMatrix2>
        void Jdiff(const Eigen::MatrixBase<StateVector1> &x0,
                   const Eigen::MatrixBase<StateVector2> &x1,
                   Eigen::MatrixBase<JMatrix1> &Jfirst, Eigen::MatrixBase<JMatrix2> &Jsecond,
                   const Jcomponent firstsecond = both) const
        {
            if (firstsecond == first || firstsecond == both)
            {
                Jfirst.setZero();
            }
            if (firstsecond == second || firstsecond == both)
            {
                Jsecond.setZero();
            }

            if (firstsecond == first)
            {
                pinocchio::dDifference(model_, head(x0, get_nq_dim()), head(x1, get_nq_dim()),
                                       topLeftCorner(Jfirst, get_nv_dim(), get_nv_dim()), pinocchio::ARG0);
                bottomRightCorner(Jfirst, get_nv_dim(), get_nv_dim()).diagonal().array() = VarScalar(-1.);
            }
            else if (firstsecond == second)
            {
                pinocchio::dDifference(model_, head(x0, get_nq_dim()), head(x1, get_nq_dim()),
                                       topLeftCorner(Jsecond, get_nv_dim(), get_nv_dim()), pinocchio::ARG1);
                bottomRightCorner(Jsecond, get_nv_dim(), get_nv_dim()).diagonal().array() = VarScalar(1.);
            }
            else
            { // computing both
                pinocchio::dDifference(model_, head(x0, get_nq_dim()), head(x1, get_nq_dim()),
                                       topLeftCorner(Jfirst, get_nv_dim(), get_nv_dim()), pinocchio::ARG0);
                pinocchio::dDifference(model_, head(x0, get_nq_dim()), head(x1, get_nq_dim()),
                                       topLeftCorner(Jsecond, get_nv_dim(), get_nv_dim()), pinocchio::ARG1);
                bottomRightCorner(Jfirst, get_nv_dim(), get_nv_dim()).diagonal().array() = VarScalar(-1.);
                bottomRightCorner(Jsecond, get_nv_dim(), get_nv_dim()).diagonal().array() = VarScalar(1.);
            }
        }

        template <typename StateVector, typename StateTangentVector, typename JMatrix1, typename JMatrix2>
        void Jintegrate(const Eigen::MatrixBase<StateVector> &x,
                        const Eigen::MatrixBase<StateTangentVector> &dx,
                        Eigen::MatrixBase<JMatrix1> &Jfirst,
                        Eigen::MatrixBase<JMatrix2> &Jsecond,
                        const Jcomponent firstsecond = both,
                        const AssignmentOp op = setto) const
        {
            if (op == setto)
            {
                if (firstsecond == first || firstsecond == both)
                    Jfirst.setZero();
                if (firstsecond == second || firstsecond == both)
                    Jsecond.setZero();
            }

            if (firstsecond == first || firstsecond == both)
            {
                switch (op)
                {
                case setto:
                    pinocchio::dIntegrate(model_, head(x, get_nq_dim()), head(dx, get_nv_dim()),
                                          topLeftCorner(Jfirst, get_nv_dim(), get_nv_dim()), pinocchio::ARG0,
                                          pinocchio::SETTO);
                    bottomRightCorner(Jfirst, get_nv_dim(), get_nv_dim()).diagonal().array() = VarScalar(1.);
                    break;
                case addto:
                    pinocchio::dIntegrate(model_, head(x, get_nq_dim()), head(dx, get_nv_dim()),
                                          topLeftCorner(Jfirst, get_nv_dim(), get_nv_dim()), pinocchio::ARG0,
                                          pinocchio::ADDTO);
                    bottomRightCorner(Jfirst, get_nv_dim(), get_nv_dim()).diagonal().array() += VarScalar(1.);
                    break;
                case rmfrom:
                    pinocchio::dIntegrate(model_, head(x, get_nq_dim()), head(dx, get_nv_dim()),
                                          topLeftCorner(Jfirst, get_nv_dim(), get_nv_dim()), pinocchio::ARG0,
                                          pinocchio::RMTO);
                    bottomRightCorner(Jfirst, get_nv_dim(), get_nv_dim()).diagonal().array() -= VarScalar(1.);
                    break;
                default:
                    break;
                }
            }
            if (firstsecond == second || firstsecond == both)
            {
                switch (op)
                {
                case setto:
                    pinocchio::dIntegrate(model_, head(x, get_nq_dim()), head(dx, get_nv_dim()),
                                          topLeftCorner(Jsecond, get_nv_dim(), get_nv_dim()), pinocchio::ARG1,
                                          pinocchio::SETTO);
                    bottomRightCorner(Jsecond, get_nv_dim(), get_nv_dim()).diagonal().array() = VarScalar(1.);
                    break;
                case addto:
                    pinocchio::dIntegrate(model_, head(x, get_nq_dim()), head(dx, get_nv_dim()),
                                          topLeftCorner(Jsecond, get_nv_dim(), get_nv_dim()), pinocchio::ARG1,
                                          pinocchio::ADDTO);
                    bottomRightCorner(Jsecond, get_nv_dim(), get_nv_dim()).diagonal().array() += VarScalar(1.);
                    break;
                case rmfrom:
                    pinocchio::dIntegrate(model_, head(x, get_nq_dim()), head(dx, get_nv_dim()),
                                          topLeftCorner(Jsecond, get_nv_dim(), get_nv_dim()), pinocchio::ARG1,
                                          pinocchio::RMTO);
                    bottomRightCorner(Jsecond, get_nv_dim(), get_nv_dim()).diagonal().array() -= VarScalar(1.);
                    break;
                default:
                    break;
                }
            }
        }

        template <typename StateVector, typename StateTangentVector, typename JMatrix>
        void JintegrateTransport(const Eigen::MatrixBase<StateVector> &x,
                                 const Eigen::MatrixBase<StateTangentVector> &dx,
                                 Eigen::MatrixBase<JMatrix> &Jin,
                                 const Jcomponent firstsecond) const
        {
            switch (firstsecond)
            {
            case first:
                pinocchio::dIntegrateTransport(model_, head(x, get_nq_dim()),
                                               head(dx, get_nv_dim()), topRows(Jin, get_nv_dim()),
                                               pinocchio::ARG0);
                break;
            case second:
                pinocchio::dIntegrateTransport(model_, head(x, get_nq_dim()),
                                               head(dx, get_nv_dim()), topRows(Jin, get_nv_dim()),
                                               pinocchio::ARG1);
                break;
            default:
                break;
            }
        }

        const RobotModel_t &get_robot() const
        {
            return model_;
        }

        const VectorNx_t &get_x0() const
        {
            return x0_;
        }

        const VectorNx_t &get_lb() const
        {
            return lb_;
        }

        const VectorNx_t &get_ub() const
        {
            return ub_;
        }

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

        using Base::get_rs;

        using Base::get_nqb;
        using Base::get_nqb_dim;

        using Base::get_nqj;
        using Base::get_nqj_dim;

        using Base::get_nq;
        using Base::get_nq_dim;

        using Base::get_nvb;
        using Base::get_nvb_dim;

        using Base::get_nvj;
        using Base::get_nvj_dim;

        using Base::get_nv;
        using Base::get_nv_dim;

        using Base::get_nrotors;
        using Base::get_nrotors_dim;

        using Base::get_nx;
        using Base::get_nx_dim;

        using Base::get_ndx;
        using Base::get_ndx_dim;

        using Base::get_nua;
        using Base::get_nua_dim;

        /**
         * @brief Display multibody state information to output stream
         */
        void disp(std::ostream &os) const
        {
            os << "StateMultibody{\n";
            os << "  RobotSpec: " << get_rs() << "\n";
            os << "  Pinocchio Model:\n";
            os << "    Model name: " << model_.name << "\n";
            os << "    Number of joints: " << model_.njoints << "\n";
            os << "    Number of bodies: " << model_.nbodies << "\n";
            os << "    Number of frames: " << model_.nframes << "\n";
            os << "    Has floating base: " << (get_nqb() > 0 ? "Yes" : "No") << "\n";
            os << "  State Bounds:\n";
            os << "    Lower bounds: " << lb_.transpose() << "\n";
            os << "    Upper bounds: " << ub_.transpose() << "\n";
            os << "    Default state: " << x0_.transpose() << "\n";
            os << "}";
        }

    protected:
        void initialize()
        {
            // Now that all of the dynamic dimensions have been updated in the constructor,
            // we can initialize the member variables
            x0_ = VectorNx_t::Zero(get_nx());
            head(x0_, get_nq_dim()) = pinocchio::neutral(model_);
            lb_ = -VectorNx_t::Constant(get_nx(), std::numeric_limits<NumScalar>::infinity());
            ub_ = VectorNx_t::Constant(get_nx(), std::numeric_limits<NumScalar>::infinity());

            head(lb_, get_nqb_dim()) = -VectorNqb_t::Constant(get_nqb(), std::numeric_limits<NumScalar>::max());
            head(ub_, get_nqb_dim()) = VectorNqb_t::Constant(get_nqb(), std::numeric_limits<NumScalar>::max());

            segment(lb_, get_nqb(), get_nqj_dim()) = tail(model_.lowerPositionLimit, get_nqj_dim());
            segment(ub_, get_nqb(), get_nqj_dim()) = tail(model_.upperPositionLimit, get_nqj_dim());

            segment(lb_, get_nq(), get_nv_dim()) = -model_.velocityLimit;
            segment(ub_, get_nq(), get_nv_dim()) = model_.velocityLimit;
        }

        RobotModel_t model_;
        VectorNx_t x0_;
        VectorNx_t lb_;
        VectorNx_t ub_;

    }; // class StateMultibodyTpl

} // namespace galileo

#endif // __galileo_core_states_state_multibody_hpp__
