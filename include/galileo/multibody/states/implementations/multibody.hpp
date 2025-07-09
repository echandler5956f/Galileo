#ifndef __galileo_multibody_states_multibody_hpp__
#define __galileo_multibody_states_multibody_hpp__

#include "galileo/core/states/state-base.hpp"
#include "galileo/multibody/robot-spec.hpp"

#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/multibody/model.hpp>

namespace galileo
{

    template <typename RobotSpec>
    class StateMultibodyTpl : public StateBase<StateMultibodyTpl<RobotSpec>, RobotSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(RS);

        StateMultibodyTpl(
            RS &rs,
            RobotModel_t *model)
            : StateBase<StateMultibodyTpl<RS>, RS>(rs),
              model_(model),
              x0_(VectorNx_t::Zero(get_nx())),
              lb_(VectorNx_t::Zero(get_nx())),
              ub_(VectorNx_t::Zero(get_nx()))
                  requires(RS::DimNQb_t::IsFixed && RS::DimNQj_t::IsFixed && RS::DimNVb_t::IsFixed && RS::DimNVj_t::IsFixed)
        {
            GALILEO_ASSERT(IsValidRobotSpec(rs), "StateMultibodyTpl: Invalid robot spec");
            initialization();
        }

        StateMultibodyTpl(
            RS &rs,
            RobotModel_t *model)
            : StateBase<StateMultibodyTpl<RS>, RS>(rs),
              model_(model)
        {
            if constexpr (RS::DimNQ_t::IsDynamic)
            {
                get_rs().NQ_dim.set_value(model->nq);
            }
            if constexpr (RS::DimNQb_t::IsDynamic)
            {
                const std::size_t nqb =
                    model->existJointName("root_joint")
                        ? model->joints[model->getJointId("root_joint")].nq()
                        : 0;
                get_rs().NQb_dim.set_value(nqb)
            }
            if constexpr (RS::DimNQj_t::IsDynamic)
            {
                const std::size_t nqj = get_nq() - get_nqb();
                get_rs().NQj_dim.set_value(nqj);
            }

            if constexpr (RS::DimNV_t::IsDynamic)
            {
                get_rs().NV_dim.set_value(model->nv);
            }
            if constexpr (RS::DimNVb_t::IsDynamic)
            {
                const std::size_t nvb =
                    model->existJointName("root_joint")
                        ? model->joints[model->getJointId("root_joint")].nv()
                        : 0;
                get_rs().NVb_dim.set_value(nvb);
            }
            if constexpr (RS::DimNVj_t::IsDynamic)
            {
                const std::size_t nvj = get_nv() - get_nvb();
                get_rs().NVj_dim.set_value(nvj);
            }

            // Now all of the dynamic dimensions are updated, so we can initialize the member variables
            x0_(VectorNx_t::Zero(get_nx()));
            lb_(VectorNx_t::Zero(get_nx()));
            ub_(VectorNx_t::Zero(get_nx()));

            GALILEO_ASSERT(IsValidRobotSpec(rs), "StateMultibodyTpl: Invalid robot spec");
            initialization();
        }

        VectorNx_t zero() const
        {
            return x0_;
        }

        VectorNx_t rand() const
        {
            VectorNx_t xrand = VectorNx_t::Random(NX);
            head(xrand, NQDim()) = pinocchio::randomConfiguration(*model_);
            return xrand;
        }

        template <typename StateVector1, typename StateVector2, typename StateTangentVector>
        void diff(const Eigen::MatrixBase<StateVector1> &x0,
                  const Eigen::MatrixBase<StateVector2> &x1,
                  Eigen::MatrixBase<StateTangentVector> &dxout) const
        {
            pinocchio::difference(*model_, head(x0, NQDim()), head(x1, NQDim()),
                                  head(dxout, NVDim()));
            tail(dxout, NVDim()) = tail(x1, NVDim()) - tail(x0, NVDim());
        }

        template <typename StateVector1, typename StateTangentVector, typename StateVector2>
        void integrate(const Eigen::MatrixBase<StateVector1> &x,
                       const Eigen::MatrixBase<StateTangentVector> &dx,
                       Eigen::MatrixBase<StateVector2> &xout) const
        {
            pinocchio::integrate(*model_, head(x, NQDim()), head(dx, NVDim()),
                                 head(xout, NQDim()));
            tail(xout, NVDim()) = tail(x, NVDim()) + tail(dx, NVDim());
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
                pinocchio::dDifference(*model_, head(x0, NQDim()), head(x1, NQDim()),
                                       topLeftCorner(Jfirst, NVDim(), NVDim()), pinocchio::ARG0);
                bottomRightCorner(Jfirst, NVDim(), NVDim()).diagonal().array() = VarScalar(-1.);
            }
            else if (firstsecond == second)
            {
                pinocchio::dDifference(*model_, head(x0, NQDim()), head(x1, NQDim()),
                                       topLeftCorner(Jsecond, NVDim(), NVDim()), pinocchio::ARG1);
                bottomRightCorner(Jsecond, NVDim(), NVDim()).diagonal().array() = VarScalar(1.);
            }
            else
            { // computing both
                pinocchio::dDifference(*model_, head(x0, NQDim()), head(x1, NQDim()),
                                       topLeftCorner(Jfirst, NVDim(), NVDim()), pinocchio::ARG0);
                pinocchio::dDifference(*model_, head(x0, NQDim()), head(x1, NQDim()),
                                       topLeftCorner(Jsecond, NVDim(), NVDim()), pinocchio::ARG1);
                bottomRightCorner(Jfirst, NVDim(), NVDim()).diagonal().array() = VarScalar(-1.);
                bottomRightCorner(Jsecond, NVDim(), NVDim()).diagonal().array() = VarScalar(1.);
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
                    pinocchio::dIntegrate(*model_, head(x, NQDim()), head(dx, NVDim()),
                                          topLeftCorner(Jfirst, NVDim(), NVDim()), pinocchio::ARG0,
                                          pinocchio::SETTO);
                    bottomRightCorner(Jfirst, NVDim(), NVDim()).diagonal().array() = VarScalar(1.);
                    break;
                case addto:
                    pinocchio::dIntegrate(*model_, head(x, NQDim()), head(dx, NVDim()),
                                          topLeftCorner(Jfirst, NVDim(), NVDim()), pinocchio::ARG0,
                                          pinocchio::ADDTO);
                    bottomRightCorner(Jfirst, NVDim(), NVDim()).diagonal().array() += VarScalar(1.);
                    break;
                case rmfrom:
                    pinocchio::dIntegrate(*model_, head(x, NQDim()), head(dx, NVDim()),
                                          topLeftCorner(Jfirst, NVDim(), NVDim()), pinocchio::ARG0,
                                          pinocchio::RMTO);
                    bottomRightCorner(Jfirst, NVDim(), NVDim()).diagonal().array() -= VarScalar(1.);
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
                    pinocchio::dIntegrate(*model_, head(x, NQDim()), head(dx, NVDim()),
                                          topLeftCorner(Jsecond, NVDim(), NVDim()), pinocchio::ARG1,
                                          pinocchio::SETTO);
                    bottomRightCorner(Jsecond, NVDim(), NVDim()).diagonal().array() = VarScalar(1.);
                    break;
                case addto:
                    pinocchio::dIntegrate(*model_, head(x, NQDim()), head(dx, NVDim()),
                                          topLeftCorner(Jsecond, NVDim(), NVDim()), pinocchio::ARG1,
                                          pinocchio::ADDTO);
                    bottomRightCorner(Jsecond, NVDim(), NVDim()).diagonal().array() += VarScalar(1.);
                    break;
                case rmfrom:
                    pinocchio::dIntegrate(*model_, head(x, NQDim()), head(dx, NVDim()),
                                          topLeftCorner(Jsecond, NVDim(), NVDim()), pinocchio::ARG1,
                                          pinocchio::RMTO);
                    bottomRightCorner(Jsecond, NVDim(), NVDim()).diagonal().array() -= VarScalar(1.);
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
                pinocchio::dIntegrateTransport(*model_, head(x, NQDim()),
                                               head(dx, NVDim()), topRows(Jin, NVDim()),
                                               pinocchio::ARG0);
                break;
            case second:
                pinocchio::dIntegrateTransport(*model_, head(x, NQDim()),
                                               head(dx, NVDim()), topRows(Jin, NVDim()),
                                               pinocchio::ARG1);
                break;
            default:
                break;
            }
        }

        const RobotModel_t *get_robot() const
        {
            return model_;
        }

        const VectorNx_t &get_x0() const
        {
            return x0_;
        }

        const VectorNx_t &get_lb_impl() const
        {
            return lb_;
        }

        const VectorNx_t &get_ub_impl() const
        {
            return ub_;
        }

        template <typename StateVector>
        void set_lb_impl(const Eigen::MatrixBase<StateVector> &lb)
        {
            lb_ = lb.derived();
        }

        template <typename StateVector>
        void set_ub_impl(const Eigen::MatrixBase<StateVector> &ub)
        {
            ub_ = ub.derived();
        }

        using Base = StateBase<StateEuclideanTpl<RS>, RS>;

        using Base::diff_dx;
        using Base::integrate_x;
        using Base::Jdiff_Js;
        using Base::Jintegrate_Js;

        using Base::get_rs;

        using Base::get_nqb;
        using Base::NQbDim;

        using Base::get_nqj;
        using Base::NQjDim;

        using Base::get_nq;
        using Base::NQDim;

        using Base::get_nvb;
        using Base::NVbDim;

        using Base::get_nvj;
        using Base::NVjDim;

        using Base::get_nv;
        using Base::NVDim;

        using Base::get_nrotors;
        using Base::NRotorsDim;

        using Base::get_nx;
        using Base::NXDim;

        using Base::get_ndx;
        using Base::NDXDim;

        using Base::get_nua;
        using Base::NUaDim;

    protected:
        void initialization()
        {
            head(x0_, NQDim()) = pinocchio::neutral(*model_);
            head(lb_, NQbDim()) = -std::numeric_limits<NumScalar>::infinity() * VectorNqb_t::Ones(get_nqb());

            head(ub_, NQbDim()) = std::numeric_limits<NumScalar>::infinity() * VectorNqb_t::Ones(get_nqb());
            segment(lb_, NQbDim(), NQjDim()) = tail(model_->lowerPositionLimit, NQjDim());
            segment(ub_, NQbDim(), NQjDim()) = tail(model_->upperPositionLimit, NQjDim());
            tail(lb_, NVDim()) = -model_->velocityLimit;
            tail(ub_, NVDim()) = model_->velocityLimit;
        }

        RobotModel_t *model_;
        VectorNx_t x0_;
        VectorNx_t lb_;
        VectorNx_t ub_;

    }; // class StateMultibodyTpl

} // namespace galileo

#endif // __galileo_multibody_states_multibody_hpp__
