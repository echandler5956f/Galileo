#ifndef __galileo_multibody_core_residuals_residual_com_position_hpp__
#define __galileo_multibody_core_residuals_residual_com_position_hpp__

#include <pinocchio/algorithm/center-of-mass-derivatives.hpp>
#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/centroidal-derivatives.hpp>
#include <pinocchio/algorithm/centroidal.hpp>

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/domains/multibody/core/residuals/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ResidualCoMPositionTpl;

    template <typename PhaseSpec>
    struct traits<ResidualCoMPositionTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualCoMPositionTpl<PS>;
        using Model_t = ResidualModelCoMPositionTpl<PS>;
        using Data_t = ResidualDataCoMPositionTpl<PS>;

        using DimNR_t = DimensionTpl<3>;
        static constexpr int NR = DimNR_t::Value;

        static constexpr bool QDependent = true;
        static constexpr bool VDependent = false;
        static constexpr bool UDependent = false;

        using R_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>>;
        using Rx_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>>;
        using Ru_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>>;
        using Arr_Rx_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>>;
        using Arr_Ru_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataCoMPositionTpl<PhaseSpec>>
    {
        using Meta_t = ResidualCoMPositionTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelCoMPositionTpl<PhaseSpec>>
    {
        using Meta_t = ResidualCoMPositionTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct ResidualDataCoMPositionTpl : public ResidualDataBase<ResidualDataCoMPositionTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ResidualCoMPositionTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualDataBase<ResidualDataCoMPositionTpl<PS>, PS>;

        using RobotData_t = typename PS::RobotData_t;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        template <typename DataCollector>
        ResidualDataCoMPositionTpl(const Model_t &model, MemoryArena &arena, DataCollector *const collector)
            : robot(collector->robot),
              R(arena, model.get_nr(), 1),
              Rx(arena, model.get_nr(), model.get_ps().get_ndx()),
              Ru(arena, model.get_nr(), model.get_ps().get_nu()),
              Arr_Rx(arena, model.get_nr(), model.get_ps().get_ndx()),
              Arr_Ru(arena, model.get_nr(), model.get_ps().get_nu())
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

    }; // class ResidualDataCoMPositionTpl

    template <typename PhaseSpec>
    class ResidualModelCoMPositionTpl : public ResidualModelBase<ResidualModelCoMPositionTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ResidualCoMPositionTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelCoMPositionTpl<PS>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        using State_t = typename PS::State_t;
        using Vector3_t = typename PS::Vector3_t;
        using VectorNu_t = typename PS::VectorNu_t;

        template <typename Vector3Type>
        ResidualModelCoMPositionTpl(const PS &ps, const State_t &state, const Eigen::MatrixBase<Vector3Type> &c_ref)
            : Base(ps, state, DimNR_t(3)), c_ref_(c_ref)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.R = data.robot->com[0] - c_ref_;
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
            leftCols(data.Rx, get_ps().get_nv_dim()) = data.robot->Jcom;
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            calcDiff(data, x, VectorNu_t::Zero(get_ps().get_nu()));
        }

        template <typename DataCollector>
        Data_t createData(MemoryArena &arena, DataCollector *const collector) const
        {
            return Data_t(*this, arena, collector);
        }

        using Base::get_ps;
        using Base::get_state;
        using Base::get_nr;
        using Base::get_nr_dim;
        using Base::get_q_dependent;
        using Base::get_u_dependent;
        using Base::get_v_dependent;

    protected:
        Vector3_t c_ref_;

    }; // class ResidualModelCoMPositionTpl

} // namespace galileo

#endif // __galileo_multibody_core_residuals_residual_com_position_hpp__
