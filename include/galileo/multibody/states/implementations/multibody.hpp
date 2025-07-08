#ifndef __galileo_multibody_states_multibody_hpp__
#define __galileo_multibody_states_multibody_hpp__

#include "galileo/core/states/state-base.hpp"

#include <pinocchio/multibody/model.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>

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
            RobotModel_t *model)
            : model_(model),
              x0_(VectorNx_t::Zero(NX))
        {
            x0_.head(NQ) = pinocchio::neutral(*model_);

            lb_.head(NQb) =
                -std::numeric_limits<NumScalar>::infinity() * VectorNqb_t::Ones(NQb);
            ub_.head(NQb) = std::numeric_limits<NumScalar>::infinity() * VectorNqb_t::Ones(NQb);
            lb_.segment(NQb, NQ - NQb) = model_->lowerPositionLimit.tail(NQ - NQb);
            ub_.segment(NQb, NQ - NQb) = model_->upperPositionLimit.tail(NQ - NQb);
            lb_.tail(NV) = -model_->velocityLimit;
            ub_.tail(NV) = model_->velocityLimit;
        }

        StateMultibodyTpl()
            : x0_(VectorNx_t::Zero(NX)) {}

        ~StateMultibodyTpl() {}

        VectorNx_t zero() const
        {
            return x0_;
        }

        VectorNx_t rand() const
        {
            VectorNx_t xrand = VectorNx_t::Random(NX);
            xrand.head(NQ) = pinocchio::randomConfiguration(*model_);
            return xrand;
        }

        template <typename StateVector1, typename StateVector2, typename StateTangentVector>
        void diff(const Eigen::MatrixBase<StateVector1> &x0,
                  const Eigen::MatrixBase<StateVector2> &x1,
                  Eigen::MatrixBase<StateTangentVector> &dxout) const
        {
            pinocchio::difference(*model_, x0.head(NQ), x1.head(NQ),
                                  dxout.head(NV));
            dxout.tail(NV) = x1.tail(NV) - x0.tail(NV);
        }

        template <typename StateVector1, typename StateTangentVector, typename StateVector2>
        void integrate(const Eigen::MatrixBase<StateVector1> &x,
                       const Eigen::MatrixBase<StateTangentVector> &dx,
                       Eigen::MatrixBase<StateVector2> &xout) const
        {
            pinocchio::integrate(*model_, x.head(NQ), dx.head(NV),
                                 xout.head(NQ));
            xout.tail(NV) = x.tail(NV) + dx.tail(NV);
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
                pinocchio::dDifference(*model_, x0.head(NQ), x1.head(NQ),
                                       Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0);
                Jfirst.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)-1;
            }
            else if (firstsecond == second)
            {
                pinocchio::dDifference(*model_, x0.head(NQ), x1.head(NQ),
                                       Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1);
                Jsecond.bottomRightCorner(NV, NV).diagonal().array() = VarScalar(1);
            }
            else
            { // computing both
                pinocchio::dDifference(*model_, x0.head(NQ), x1.head(NQ),
                                       Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0);
                pinocchio::dDifference(*model_, x0.head(NQ), x1.head(NQ),
                                       Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1);
                Jfirst.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)-1;
                Jsecond.bottomRightCorner(NV, NV).diagonal().array() = (VarScalar)1;
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
                    pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                          Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0,
                                          pinocchio::SETTO);
                    Jfirst.bottomRightCorner(NV, NV).diagonal().array() = VarScalar(1);
                    break;
                case addto:
                    pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                          Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0,
                                          pinocchio::ADDTO);
                    Jfirst.bottomRightCorner(NV, NV).diagonal().array() += VarScalar(1);
                    break;
                case rmfrom:
                    pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                          Jfirst.topLeftCorner(NV, NV), pinocchio::ARG0,
                                          pinocchio::RMTO);
                    Jfirst.bottomRightCorner(NV, NV).diagonal().array() -= VarScalar(1);
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
                    pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                          Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1,
                                          pinocchio::SETTO);
                    Jsecond.bottomRightCorner(NV, NV).diagonal().array() = VarScalar(1);
                    break;
                case addto:
                    pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                          Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1,
                                          pinocchio::ADDTO);
                    Jsecond.bottomRightCorner(NV, NV).diagonal().array() += VarScalar(1);
                    break;
                case rmfrom:
                    pinocchio::dIntegrate(*model_, x.head(NQ), dx.head(NV),
                                          Jsecond.topLeftCorner(NV, NV), pinocchio::ARG1,
                                          pinocchio::RMTO);
                    Jsecond.bottomRightCorner(NV, NV).diagonal().array() -= VarScalar(1);
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
                pinocchio::dIntegrateTransport(*model_, x.head(NQ),
                                               dx.head(NV), Jin.topRows(NV),
                                               pinocchio::ARG0);
                break;
            case second:
                pinocchio::dIntegrateTransport(*model_, x.head(NQ),
                                               dx.head(NV), Jin.topRows(NV),
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

        /**
         * @brief Return the state lower bound
         */
        const VectorNx_t &get_lb() const
        {
            return lb_;
        }

        /**
         * @brief Return the state upper bound
         */
        const VectorNx_t &get_ub() const
        {
            return ub_;
        }

        /**
         * @brief Modify the state lower bound
         */
        template <typename StateVector>
        void set_lb(const Eigen::MatrixBase<StateVector> &lb)
        {
            lb_ = lb.derived();
        }

        /**
         * @brief Modify the state upper bound
         */
        template <typename StateVector>
        void set_ub(const Eigen::MatrixBase<StateVector> &ub)
        {
            ub_ = ub.derived();
        }

    protected:
        RobotModel_t *model_;
        VectorNx_t x0_;
        VectorNx_t lb_;
        VectorNx_t ub_;

    }; // class StateMultibodyTpl

} // namespace galileo

#endif // __galileo_multibody_states_multibody_hpp__
