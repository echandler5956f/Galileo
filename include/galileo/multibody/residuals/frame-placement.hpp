#ifndef __galileo_multibody_residuals_frame_placement_hpp__
#define __galileo_multibody_residuals_frame_placement_hpp__

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/frames-derivatives.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>

#include "galileo/core/residuals/residual-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ResidualFramePlacementTpl;

    template <typename PhaseSpec>
    struct traits<ResidualFramePlacementTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualFramePlacementTpl<PS>;
        using Model_t = ResidualModelFramePlacementTpl<PS>;
        using Data_t = ResidualDataFramePlacementTpl<PS>;

        static constexpr int NR = 6;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataFramePlacementTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualFramePlacementTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelFramePlacementTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualFramePlacementTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct ResidualDataFramePlacementTpl : public ResidualDataBase<ResidualDataFramePlacementTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ResidualFramePlacementTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        ResidualDataFramePlacementTpl() : R(R_t::Zero()), Rx(Rx_t::Zero()), Ru(Ru_t::Zero()), Arr_Rx(Arr_Rx_t::Zero()), Arr_Ru(Arr_Ru_t::Zero())
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

        SE3_t rMf;
        Matrix6_t rJf;
        Matrix6Nv_t fJf;

    }; // class ResidualDataFramePlacementTpl

    template <typename PhaseSpec>
    class ResidualModelFramePlacementTpl : public ResidualModelBase<ResidualModelFramePlacementTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ResidualFramePlacementTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using RobotModel_t = typename PS::RobotModel_t;
        using FrameIndex_t = pinocchio::FrameIndex;
        using SE3_t = pinocchio::SE3Tpl<typename PS::VarScalar>;

        ResidualModelFramePlacementTpl(RobotModel_t *robot_model,
                                       const FrameIndex_t frame_id,
                                       const SE3_t &p_ref,
                                       const int nu)
            : robot_model_(robot_model), frame_id_(frame_id), p_ref_(p_ref), oMf_inv_(p_ref.inverse()), nu_(nu)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            pinocchio::updateFramePlacement(*robot_model_, *data.robot, frame_id_);
            data.rMf = oMf_inv_ * data.robot->oMf[frame_id_];
            data.R = pinocchio::log6(data.rMf).toVector();
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            pinocchio::Jlog6(data.rMf, data.rJf);
            pinocchio::getFrameJacobian(
                *robot_model_,
                *data.robot,
                frame_id_,
                pinocchio::ReferenceFrame::LOCAL,
                data.fJf);
            data.Rx.leftCols(PS::NV).noalias() = data.rJf * data.fJf;
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
        SE3_t p_ref_;
        SE3_t oMf_inv_;
        int nu_;

    }; // class ResidualModelFramePlacementTpl

} // namespace galileo

#endif // __galileo_multibody_residuals_frame_placement_hpp__
