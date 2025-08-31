#ifndef __galileo_core_residuals_residual_control_hpp__
#define __galileo_core_residuals_residual_control_hpp__

#include "galileo/core/residuals/residual-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ResidualControlTpl;

    template <typename PhaseSpec>
    struct traits<ResidualControlTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualControlTpl<PS>;
        using Model_t = ResidualModelControlTpl<PS>;
        using Data_t = ResidualDataControlTpl<PS>;

        using DimNR_t = typename PS::DimNU_t;
        static constexpr int NR = DimNR_t::Value;

        static constexpr bool QDependent = false;
        static constexpr bool VDependent = false;
        static constexpr bool UDependent = true;

        using R_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>>;
        using Rx_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>>;
        using Ru_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>>;
        using Arr_Rx_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>>;
        using Arr_Ru_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataControlTpl<PhaseSpec>>
    {
        using Meta_t = ResidualControlTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelControlTpl<PhaseSpec>>
    {
        using Meta_t = ResidualControlTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct ResidualDataControlTpl : public ResidualDataBase<ResidualDataControlTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ResidualControlTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualDataBase<ResidualDataControlTpl<PS>, PS>;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        template <typename DataCollector>
        ResidualDataControlTpl(const Model_t &model, MemoryArena &arena, DataCollector *const collector)
            : R(arena, model.get_nr(), 1),
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
            Ru.diagonal().setOnes();
        }

        R_t R;
        Rx_t Rx;
        Ru_t Ru;
        Arr_Rx_t Arr_Rx;
        Arr_Ru_t Arr_Ru;

    }; // class ResidualDataControlTpl

    template <typename PhaseSpec>
    class ResidualModelControlTpl : public ResidualModelBase<ResidualModelControlTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ResidualControlTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelControlTpl<PS>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        using State_t = typename PS::State_t;
        using VectorNu_t = typename PS::VectorNU_t;

        template <typename ControlVectorType>
        ResidualModelControlTpl(const PS &ps, const State_t &state, const Eigen::MatrixBase<ControlVectorType> &u_ref)
            : Base(ps, state, DimNR_t(ps.get_nu())), u_ref_(u_ref)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.R = u - u_ref_;
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            data.R.setZero();
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            // The Jacobian has constant values which were set in createData
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            // The Jacobian has constant values which were set in createData
        }

        template <bool UpdateU = true, typename CostDataType, typename ActivationDataType>
        void calcCostDiffImpl(CostDataType &cdata, Data_t &rdata, const ActivationDataType &adata) const
        {
            cdata.Lu = adata.Ar;
            cdata.Luu = adata.Arr;
        }

        template <typename DataCollector>
        Data_t createData(MemoryArena &arena, DataCollector *const collector) const
        {
            return Data_t(*this, arena,  collector);
        }

        using Base::get_ps;
        using Base::get_state;
        using Base::get_nr;
        using Base::get_nr_dim;
        using Base::get_q_dependent;
        using Base::get_u_dependent;
        using Base::get_v_dependent;

    protected:
        VectorNu_t u_ref_;

    }; // class ResidualModelControlTpl

} // namespace galileo

#endif // __galileo_core_residuals_residual_control_hpp__
