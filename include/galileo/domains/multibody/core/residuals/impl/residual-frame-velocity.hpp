#ifndef __galileo_multibody_core_residuals_residual_frame_velocity_hpp__
#define __galileo_multibody_core_residuals_residual_frame_velocity_hpp__

#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/frames-derivatives.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/domains/multibody/core/residuals/fwd.hpp"

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
        static constexpr int NR = DimNR_t::Value;

        static constexpr bool QDependent = true;
        static constexpr bool VDependent = true;
        static constexpr bool UDependent = false;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataFrameVelocityTpl<PhaseSpec>>
    {
        using Meta_t = ResidualFrameVelocityTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelFrameVelocityTpl<PhaseSpec>>
    {
        using Meta_t = ResidualFrameVelocityTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct ResidualDataFrameVelocityTpl : public ResidualDataBase<ResidualDataFrameVelocityTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualDataBase<ResidualDataFrameVelocityTpl<PS>, PS>;

        using RobotData_t = typename PS::RobotData_t;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        template <typename DataCollector>
        ResidualDataFrameVelocityTpl(const Model_t &model, DataCollector *const collector)
            : robot(collector->robot),
              R(model.get_nr()),
              Rx(model.get_nr(), model.get_ps().get_ndx()),
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

        std::shared_ptr<RobotData_t> robot;

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
        using PS = PhaseSpec;

        using Meta_t = ResidualFrameVelocityTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelFrameVelocityTpl<PS>, PS>;

        using State_t = typename PS::State_t;
        using VectorNu_t = typename PS::VectorNu_t;
        using FrameIndex_t = typename PS::FrameIndex_t;
        using Motion_t = typename PS::Motion_t;
        using ReferenceFrame_t = typename PS::ReferenceFrame_t;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        ResidualModelFrameVelocityTpl(const PS &ps,
                                      const State_t &state,
                                      const FrameIndex_t frame_id,
                                      const Motion_t &velocity,
                                      const ReferenceFrame_t type)
            : Base(ps, state, DimNR_t(6)), frame_id_(frame_id), vref_(velocity), type_(type)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.R =
                (pinocchio::getFrameVelocity(get_state().get_robot(), *data.robot.get(), frame_id_, type_) -
                 vref_)
                    .toVector();
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            calc(data, x, VectorNu_t::Zero(get_ps().get_nu()));
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            pinocchio::getFrameVelocityDerivatives(get_state().get_robot(),
                                                   *data.robot.get(),
                                                   frame_id_,
                                                   type_,
                                                   leftCols(data.Rx, get_ps().get_nv_dim()),
                                                   rightCols(data.Rx, get_ps().get_nv_dim()));
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            calcDiff(data, x, VectorNu_t::Zero(get_ps().get_nu()));
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return Data_t(*this, collector);
        }

        using Base::get_ps;
        using Base::get_state;
        using Base::get_nr;
        using Base::get_nr_dim;
        using Base::get_q_dependent;
        using Base::get_u_dependent;
        using Base::get_v_dependent;

    protected:
        FrameIndex_t frame_id_;
        Motion_t vref_;
        ReferenceFrame_t type_;

    }; // class ResidualModelFrameVelocityTpl

} // namespace galileo

#endif // __galileo_multibody_core_residuals_residual_frame_velocity_hpp__
