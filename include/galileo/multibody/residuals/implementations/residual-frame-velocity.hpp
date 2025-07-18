#ifndef __galileo_multibody_residuals_residual_frame_velocity_hpp__
#define __galileo_multibody_residuals_residual_frame_velocity_hpp__

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/frames-derivatives.hpp>
#include <pinocchio/algorithm/frames.hpp>
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

        using DimNR_t = DimensionTpl<6>;

        static constexpr bool QDependent = true;
        static constexpr bool VDependent = true;
        static constexpr bool UDependent = false;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNU_t::Value, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNU_t::Value, PS::Options>;
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
    struct ResidualDataFrameVelocityTpl
        : public ResidualDataBase<ResidualDataFrameVelocityTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualDataBase<ResidualDataFrameVelocityTpl<PS>, PS>;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        template <typename DataCollector>
        ResidualDataFrameVelocityTpl(const Model_t &model, DataCollector *const collector)
            : R(model.get_nr()), Rx(model.get_nr(), model.get_ps().get_ndx()),
              Ru(model.get_nr(), model.get_ps().get_nu()),
              Arr_Rx(model.get_nr(), model.get_ps().get_ndx()),
              Arr_Ru(model.get_nr(), model.get_ps().get_nu())
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
    class ResidualModelFrameVelocityTpl
        : public ResidualModelBase<ResidualModelFrameVelocityTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelFrameVelocityTpl<PS>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        using FrameIndex_t = pinocchio::FrameIndex;
        using Motion_t = pinocchio::MotionTpl<typename PS::NumScalar>;
        using ReferenceFrame_t = pinocchio::ReferenceFrame;

        ResidualModelFrameVelocityTpl(const PS &ps,
                                      const std::shared_ptr<State_t> &state,
                                      const FrameIndex_t frame_id,
                                      const Motion_t &velocity,
                                      const ReferenceFrame_t type)
            : Base(ps, DimNR_t()),
              state_(state), frame_id_(frame_id), vref_(velocity), type_(type)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.R = (pinocchio::getFrameVelocity(
                          state_->get_robot(),
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
                state_->get_robot(),
                *data.robot,
                frame_id_,
                type_,
                leftCols(data.Rx, get_ps().get_nv_dim()),
                rightCols(data.Rx, get_ps().get_nv_dim()));
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            Data_t data(*this, collector);
            data.robot = collector->robot;
            return data;
        }

        const std::shared_ptr<State_t> &get_state() const
        {
            return state_;
        }

        using Base::get_ps;

        using Base::get_nr;
        using Base::get_nr_dim;

        using Base::get_q_dependent;
        using Base::get_v_dependent;
        using Base::get_u_dependent;

    protected:
        std::shared_ptr<State_t> state_;
        FrameIndex_t frame_id_;
        Motion_t vref_;
        ReferenceFrame_t type_;

    }; // class ResidualModelFrameVelocityTpl

} // namespace galileo

#endif // __galileo_multibody_residuals_residual_frame_velocity_hpp__
