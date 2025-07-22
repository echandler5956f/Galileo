#ifndef __galileo_multibody_residuals_residual_control_hpp__
#define __galileo_multibody_residuals_residual_control_hpp__

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/multibody/residuals/fwd.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

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

        static constexpr bool QDependent = false;
        static constexpr bool VDependent = false;
        static constexpr bool UDependent = true;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNU_t::Value, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNU_t::Value, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataControlTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualControlTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelControlTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualControlTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct ResidualDataControlTpl
        : public ResidualDataBase<ResidualDataControlTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

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
        ResidualDataControlTpl(const Model_t &model, DataCollector *const collector)
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
            Ru.diagonal().fill(1.0);
        }

        R_t R;
        Rx_t Rx;
        Ru_t Ru;
        Arr_Rx_t Arr_Rx;
        Arr_Ru_t Arr_Ru;

    }; // class ResidualDataControlTpl

    template <typename PhaseSpec>
    class ResidualModelControlTpl
        : public ResidualModelBase<ResidualModelControlTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = ResidualControlTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelControlTpl<PS>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;


        template <typename ControlVectorType>
        ResidualModelControlTpl(const PS &ps,
                                const Eigen::MatrixBase<ControlVectorType> &u_ref)
            : Base(ps, DimNR_t()),
              u_ref_(u_ref)
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
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
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
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            // The Jacobian has constant values which were set in createData
        }

        template <typename CostDataType, typename ActivationDataType, bool UpdateU = true>
        void calcCostDiffImpl(CostDataType &cdata,
                              Data_t &rdata,
                              const ActivationDataType &adata) const
        {
            cdata.Lu = adata.Ar;
            cdata.Luu = adata.Arr;
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return Data_t(*this, collector);
        }

        using Base::get_ps;

        using Base::get_nr;
        using Base::get_nr_dim;

        using Base::get_q_dependent;
        using Base::get_u_dependent;
        using Base::get_v_dependent;

    protected:
        VectorNu_t u_ref_;

    }; // class ResidualModelControlTpl

} // namespace galileo

#endif // __galileo_multibody_residuals_residual_control_hpp__
