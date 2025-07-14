#ifndef __galileo_core_costs_cost_residual_hpp__
#define __galileo_core_costs_cost_residual_hpp__

#include "galileo/core/costs/cost-base.hpp"

namespace galileo
{

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        template <typename PS> class ActivationTpl>
    struct CostResidualTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        template <typename PS> class ActivationTpl>
    struct traits<CostResidualTpl<PhaseSpec, ResidualTpl, ActivationTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = CostResidualTpl<PS, ResidualTpl, ActivationTpl>;
        using Model_t = CostModelResidualTpl<PS, ResidualTpl, ActivationTpl>;
        using Data_t = CostDataResidualTpl<PS, ResidualTpl, ActivationTpl>;

        using ResidualMeta_t = typename traits<ResidualTpl<PS>>::Meta_t;
        using ResidualModel_t = typename traits<ResidualMeta_t>::Model_t;
        using ResidualData_t = typename traits<ResidualMeta_t>::Data_t;

        using ActivationMeta_t = typename traits<ActivationTpl<PS>>::Meta_t;
        using ActivationModel_t = typename traits<ActivationMeta_t>::Model_t;
        using ActivationData_t = typename traits<ActivationMeta_t>::Data_t;

        using DimNR_t = typename traits<ResidualMeta_t>::DimNR_t;

        using L_t = typename PS::VarScalar;
        using Lx_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNDX_t::Value, 1, PS::Options>;
        using Lu_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNU_t::Value, 1, PS::Options>;
        using Lxx_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNDX_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Lxu_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNDX_t::Value, PS::DimNU_t::Value, PS::Options>;
        using Luu_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNU_t::Value, PS::DimNU_t::Value, PS::Options>;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        template <typename PS> class ActivationTpl>
    struct traits<CostDataResidualTpl<PhaseSpec, ResidualTpl, ActivationTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = CostResidualTpl<PS, ResidualTpl, ActivationTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        template <typename PS> class ActivationTpl>
    struct traits<CostModelResidualTpl<PhaseSpec, ResidualTpl, ActivationTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = CostResidualTpl<PS, ResidualTpl, ActivationTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        template <typename PS> class ActivationTpl>
    struct CostDataResidualTpl
        : public CostDataBase<CostDataResidualTpl<PhaseSpec, ResidualTpl, ActivationTpl>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = CostResidualTpl<PS, ResidualTpl, ActivationTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Base = CostDataBase<CostDataResidualTpl<PhaseSpec, ResidualTpl, ActivationTpl>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        GALILEO_COST_DATA_TYPEDEF(Meta_t);

        using ResidualMeta_t = typename traits<Meta_t>::ResidualMeta_t;
        using ResidualModel_t = typename traits<Meta_t>::ResidualModel_t;
        using ResidualData_t = typename traits<Meta_t>::ResidualData_t;

        using ActivationMeta_t = typename traits<Meta_t>::ActivationMeta_t;
        using ActivationModel_t = typename traits<Meta_t>::ActivationModel_t;
        using ActivationData_t = typename traits<Meta_t>::ActivationData_t;

        DEFAULT_ACCESSOR(L_t, L);
        DEFAULT_ACCESSOR(Lx_t, Lx);
        DEFAULT_ACCESSOR(Lu_t, Lu);
        DEFAULT_ACCESSOR(Lxx_t, Lxx);
        DEFAULT_ACCESSOR(Lxu_t, Lxu);
        DEFAULT_ACCESSOR(Luu_t, Luu);

        ResidualData_t residual;
        ActivationData_t activation;
        L_t L;
        Lx_t Lx;
        Lu_t Lu;
        Lxx_t Lxx;
        Lxu_t Lxu;
        Luu_t Luu;

    }; // struct CostDataResidualTpl

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        template <typename PS> class ActivationTpl>
    class CostModelResidualTpl
        : public CostModelBase<CostModelResidualTpl<PhaseSpec, ResidualTpl, ActivationTpl>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = CostResidualTpl<PS, ResidualTpl, ActivationTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Base = CostModelBase<CostModelResidualTpl<PhaseSpec, ResidualTpl, ActivationTpl>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        using ResidualMeta_t = typename traits<Meta_t>::ResidualMeta_t;
        using ResidualModel_t = typename traits<Meta_t>::ResidualModel_t;
        using ResidualData_t = typename traits<Meta_t>::ResidualData_t;

        using ActivationMeta_t = typename traits<Meta_t>::ActivationMeta_t;
        using ActivationModel_t = typename traits<Meta_t>::ActivationModel_t;
        using ActivationData_t = typename traits<Meta_t>::ActivationData_t;

        CostModelResidualTpl(const ResidualModel_t &residual,
                             const ActivationModel_t &activation)
            : residual_(residual),
              activation_(activation)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            residual_.calc(data.residual, x.derived(), u.derived());
            activation_.calc(data.activation, x.derived(), u.derived());

            data.L = data.activation.A;
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            residual_.calc(data.residual, x.derived());
            activation_.calc(data.activation, x.derived());

            data.L = data.activation.A;
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            residual_.calcDiff(data.residual, x.derived(), u.derived());
            activation_.calcDiff(data.activation, x.derived(), u.derived());

            residual_.calcCostDiff(data, data.residual, data.activation);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            residual_.calcDiff(data.residual, x.derived());
            activation_.calcDiff(data.activation, x.derived());

            residual_.calcCostDiff(data, data.residual, data.activation, false);
        }

    protected:
        ResidualModel_t residual_;
        ActivationModel_t activation_;

    }; // class CostModelResidualTpl

} // namespace galileo

#endif // __galileo_core_costs_cost_residual_hpp__
