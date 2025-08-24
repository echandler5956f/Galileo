#ifndef __galileo_multibody_core_residuals_residual_frame_translation_hpp__
#define __galileo_multibody_core_residuals_residual_frame_translation_hpp__

#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/frames-derivatives.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/domains/multibody/core/residuals/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ResidualFrameTranslationTpl;

    template <typename PhaseSpec>
    struct traits<ResidualFrameTranslationTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualFrameTranslationTpl<PS>;
        using Model_t = ResidualModelFrameTranslationTpl<PS>;
        using Data_t = ResidualDataFrameTranslationTpl<PS>;

        using DimNR_t = DimensionTpl<3>;
        static constexpr int NR = DimNR_t::Value;

        static constexpr bool QDependent = true;
        static constexpr bool VDependent = false;
        static constexpr bool UDependent = false;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataFrameTranslationTpl<PhaseSpec>>
    {
        using Meta_t = ResidualFrameTranslationTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelFrameTranslationTpl<PhaseSpec>>
    {
        using Meta_t = ResidualFrameTranslationTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct ResidualDataFrameTranslationTpl
        : public ResidualDataBase<ResidualDataFrameTranslationTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ResidualFrameTranslationTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualDataBase<ResidualDataFrameTranslationTpl<PS>, PS>;

        using RobotData_t = typename PS::RobotData_t;
        using Matrix6Nv_t = typename PS::Matrix6Nv_t;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        template <typename DataCollector>
        ResidualDataFrameTranslationTpl(const Model_t &model, DataCollector *const collector)
            : robot(collector->robot),
              R(model.get_nr()),
              Rx(model.get_nr(), model.get_ps().get_ndx()),
              Ru(model.get_nr(), model.get_ps().get_nu()),
              Arr_Rx(model.get_nr(), model.get_ps().get_ndx()),
              Arr_Ru(model.get_nr(), model.get_ps().get_nu()),
              fJf(6, model.get_ps().get_nv())
        {
            R.setZero();
            Rx.setZero();
            Ru.setZero();
            Arr_Rx.setZero();
            Arr_Ru.setZero();

            fJf.setZero();
        }

        std::shared_ptr<RobotData_t> robot;

        R_t R;
        Rx_t Rx;
        Ru_t Ru;
        Arr_Rx_t Arr_Rx;
        Arr_Ru_t Arr_Ru;

        Matrix6Nv_t fJf;

    }; // class ResidualDataFrameTranslationTpl

    template <typename PhaseSpec>
    class ResidualModelFrameTranslationTpl
        : public ResidualModelBase<ResidualModelFrameTranslationTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ResidualFrameTranslationTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelFrameTranslationTpl<PS>, PS>;

        using State_t = typename PS::State_t;
        using Vector3_t = typename PS::Vector3_t;
        using VectorNu_t = typename PS::VectorNu_t;
        using FrameIndex_t = typename PS::FrameIndex_t;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        template <typename Vector3Type>
        ResidualModelFrameTranslationTpl(const PS &ps,
                                         const State_t &state,
                                         const FrameIndex_t frame_id,
                                         const Eigen::MatrixBase<Vector3Type> &x_ref)
            : Base(ps, state, DimNR_t(3)), frame_id_(frame_id), x_ref_(x_ref)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            pinocchio::updateFramePlacement(get_state().get_robot(), *data.robot.get(), frame_id_);
            data.R = data.robot->oMf[frame_id_].translation() - x_ref_;
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
            pinocchio::getFrameJacobian(get_state().get_robot(),
                                        *data.robot.get(),
                                        frame_id_,
                                        pinocchio::ReferenceFrame::LOCAL,
                                        data.fJf);

            leftCols(data.Rx, get_ps().get_nv_dim()).noalias() =
                data.robot->oMf[frame_id_].rotation() * topRows(data.fJf, 3);
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
        Vector3_t x_ref_;

    }; // class ResidualModelFrameTranslationTpl

} // namespace galileo

#endif // __galileo_multibody_core_residuals_residual_frame_translation_hpp__
