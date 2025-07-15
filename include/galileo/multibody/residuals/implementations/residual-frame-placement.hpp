#ifndef __galileo_multibody_residuals_residual_frame_placement_hpp__
#define __galileo_multibody_residuals_residual_frame_placement_hpp__

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/frames-derivatives.hpp>
#include <pinocchio/algorithm/frames.hpp>
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

        using DimNR_t = DimensionTpl<6>;
        static constexpr int NR = DimNR_t::Value;

        static constexpr bool QDependent = true;
        static constexpr bool VDependent = false;
        static constexpr bool UDependent = false;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::DimNDX_t::Value, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::DimNU_t::Value, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::DimNDX_t::Value, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::DimNU_t::Value, PS::Options>;
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
    struct ResidualDataFramePlacementTpl
        : public ResidualDataBase<ResidualDataFramePlacementTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ResidualFramePlacementTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualDataBase<ResidualDataFramePlacementTpl<PS>, PS>;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        ResidualDataFramePlacementTpl(const Model_t &model)
            : R(model.get_nr()), Rx(model.get_nr(), model.get_ps().ndx_dim.value()),
              Ru(model.get_nr(), model.get_ps().nu_dim.value()),
              Arr_Rx(model.get_nr(), model.get_ps().ndx_dim.value()),
              Arr_Ru(model.get_nr(), model.get_ps().nu_dim.value()),
              rJf(6, 6), fJf(6, model.get_ps().nv_dim.value())
        {
            R.setZero();
            Rx.setZero();
            Ru.setZero();
            Arr_Rx.setZero();
            Arr_Ru.setZero();

            rJf.setZero();
            fJf.setZero();
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
    class ResidualModelFramePlacementTpl
        : public ResidualModelBase<ResidualModelFramePlacementTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ResidualFramePlacementTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelFramePlacementTpl<PS>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        using FrameIndex_t = pinocchio::FrameIndex;
        using SE3_t = pinocchio::SE3Tpl<typename PS::VarScalar>;

        ResidualModelFramePlacementTpl(const PS &ps,
                                       const std::shared_ptr<State_t> &state,
                                       const FrameIndex_t frame_id,
                                       const SE3_t &p_ref)
            : Base(ps, DimNR_t()),
              state_(state), frame_id_(frame_id), p_ref_(p_ref), oMf_inv_(p_ref.inverse())
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            pinocchio::updateFramePlacement(state_->get_robot(), *data.robot, frame_id_);
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
                state_->get_robot(),
                *data.robot,
                frame_id_,
                pinocchio::ReferenceFrame::LOCAL,
                data.fJf);
            leftCols(data.Rx, get_ps().nv_dim).noalias() = data.rJf * data.fJf;
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            Data_t data(*this);
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

        using Base::get_nu;
        using Base::get_nu_dim;

        using Base::get_q_dependent;
        using Base::get_v_dependent;
        using Base::get_u_dependent;

    protected:
        std::shared_ptr<State_t> state_;
        FrameIndex_t frame_id_;
        SE3_t p_ref_;
        SE3_t oMf_inv_;

    }; // class ResidualModelFramePlacementTpl

} // namespace galileo

#endif // __galileo_multibody_residuals_residual_frame_placement_hpp__
