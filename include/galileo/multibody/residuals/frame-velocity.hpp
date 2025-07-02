#ifndef __galileo_multibody_residuals_frame_velocity_hpp__
#define __galileo_multibody_residuals_frame_velocity_hpp__

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/frames-derivatives.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>

#include "galileo/core/residuals/residual-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ResidualFrameVelocityTpl;

    template <typename PhaseSpec>
    struct traits<ResidualFrameVelocityTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = ResidualModelFrameVelocityTpl<PS>;
        using Data_t = ResidualDataFrameVelocityTpl<PS>;

        static constexpr int NR = 6;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataFrameVelocityTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelFrameVelocityTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct ResidualDataFrameVelocityTpl : public ResidualDataBase<ResidualDataFrameVelocityTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        ResidualDataFrameVelocityTpl() : R(R_t::Zero()), Rx(Rx_t::Zero()), Ru(Ru_t::Zero()), Arr_Rx(Arr_Rx_t::Zero()), Arr_Ru(Arr_Ru_t::Zero())
        {
            R.setZero();
            Rx.setZero();
            Ru.setZero();
            Arr_Rx.setZero();
            Arr_Ru.setZero();
        }

        typename PS::RobotData_t *robot;

        R_t R;
        Rx_t Rx;
        Ru_t Ru;
        Arr_Rx_t Arr_Rx;
        Arr_Ru_t Arr_Ru;

    }; // class ResidualDataFrameVelocityTpl

    template <typename PhaseSpec>
    class ResidualModelFrameVelocityTpl : public ResidualModelBase<ResidualModelFrameVelocityTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using RobotModel_t = typename PS::RobotModel_t;
        using FrameIndex_t = pinocchio::FrameIndex;
        using Motion_t = pinocchio::MotionTpl<typename PS::NumScalar>;
        using ReferenceFrame_t = pinocchio::ReferenceFrame;

        ResidualModelFrameVelocityTpl(RobotModel_t *robot_model,
                                      const FrameIndex_t frame_id,
                                      const Motion_t &velocity,
                                      const ReferenceFrame_t type,
                                      const int nu)
            : robot_model_(robot_model), frame_id_(frame_id), vref_(velocity), type_(type), nu_(nu)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.R = (pinocchio::getFrameVelocity(
                          *robot_model_,
                          *data.robot,
                          frame_id_,
                          type_) -
                      vref_)
                         .toVector();
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            pinocchio::getFrameVelocityDerivatives(
                *robot_model_,
                *data.robot,
                frame_id_,
                type_,
                data.Rx.leftCols(PS::NV),
                data.Rx.rightCols(PS::NV));
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            Data_t data;
            data.robot = collector->robot;
            return data;
        }

        bool q_dependent_impl() const
        {
            return true;
        }

        bool v_dependent_impl() const
        {
            return true;
        }

        bool u_dependent_impl() const
        {
            return false;
        }

        int nu_impl() const
        {
            return nu_;
        }

    protected:
        RobotModel_t *robot_model_;
        FrameIndex_t frame_id_;
        Motion_t vref_;
        ReferenceFrame_t type_;
        int nu_;

    }; // class ResidualModelFrameVelocityTpl

} // namespace galileo

#endif // __galileo_multibody_residuals_frame_velocity_hpp__
