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
        using Base = CostDataBase<CostDataResidualTpl<PS, ResidualTpl, ActivationTpl>, PS>;

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

        template <typename DataCollector>
        CostDataResidualTpl(const Model_t &model, DataCollector *const collector)
            : activation(model.get_activation().createData()),
              residual(model.get_residual().createData(collector)),
              L(L_t(0.)),
              Lx(model.get_ps().get_ndx()),
              Lu(model.get_ps().get_nu()),
              Lxx(model.get_ps().get_ndx(), model.get_ps().get_ndx()),
              Lxu(model.get_ps().get_ndx(), model.get_ps().get_nu()),
              Luu(model.get_ps().get_nu(), model.get_ps().get_nu())
        {
            Lx.setZero();
            Lu.setZero();
            Lxx.setZero();
            Lxu.setZero();
            Luu.setZero();
        }

        ActivationData_t activation;
        ResidualData_t residual;
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
        using Base = CostModelBase<CostModelResidualTpl<PS, ResidualTpl, ActivationTpl>, PS>;

        using ResidualMeta_t = typename traits<Meta_t>::ResidualMeta_t;
        using ResidualModel_t = typename traits<Meta_t>::ResidualModel_t;
        using ResidualData_t = typename traits<Meta_t>::ResidualData_t;

        using ActivationMeta_t = typename traits<Meta_t>::ActivationMeta_t;
        using ActivationModel_t = typename traits<Meta_t>::ActivationModel_t;
        using ActivationData_t = typename traits<Meta_t>::ActivationData_t;

        CostModelResidualTpl(const PS &ps, const ResidualModel_t &residual,
                             const ActivationModel_t &activation)
            : Base(ps),
              residual_(residual),
              activation_(activation)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            std::cout << "residual calc" << std::endl;
            residual_.calc(data.residual, x.derived(), u.derived());
            std::cout << "residual calc done" << std::endl;
            activation_.calc(data.activation, data.residual.R);
            std::cout << "activation calc done" << std::endl;
            data.L = data.activation.A;
            std::cout << "data.L: " << data.L << std::endl;
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            if constexpr (!residual_.get_q_dependent() && !residual_.get_v_dependent())
            {
                data.activation.A = 0.;
                data.L = 0.;
                return;
            }
            else
            {
                residual_.calc(data.residual, x.derived());
                activation_.calc(data.activation, data.residual.R);

                data.L = data.activation.A;
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            residual_.calcDiff(data.residual, x.derived(), u.derived());
            activation_.calcDiff(data.activation, data.residual.R);

            residual_.calcCostDiff<true>(data, data.residual, data.activation);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            if constexpr (!residual_.get_q_dependent() && !residual_.get_v_dependent())
            {
                data.Lx.setZero();
                data.Lxx.setZero();
                return;
            }
            else
            {
                residual_.calcDiff(data.residual, x.derived());
                activation_.calcDiff(data.activation, data.residual.R);

                residual_.calcCostDiff<false>(data, data.residual, data.activation);
            }
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return Data_t(*this, collector);
        }

        const ResidualModel_t &get_residual() const
        {
            return residual_;
        }

        const ActivationModel_t &get_activation() const
        {
            return activation_;
        }

        using Base::get_ps;

    protected:
        ResidualModel_t residual_;
        ActivationModel_t activation_;

    }; // class CostModelResidualTpl

} // namespace galileo

#endif // __galileo_core_costs_cost_residual_hpp__
