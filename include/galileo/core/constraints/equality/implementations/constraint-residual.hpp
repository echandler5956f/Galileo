#ifndef __galileo_core_constraints_equality_constraint_residual_hpp__
#define __galileo_core_constraints_equality_constraint_residual_hpp__

#include "galileo/core/constraints/equality/constraint-base.hpp"

namespace galileo
{

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl>
    struct ConstraintResidualTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl>
    struct traits<ConstraintResidualTpl<PhaseSpec, ResidualTpl>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = ConstraintResidualTpl<PS, ResidualTpl>;
        using Model_t = ConstraintModelResidualTpl<PS, ResidualTpl>;
        using Data_t = ConstraintDataResidualTpl<PS, ResidualTpl>;

        using ResidualMeta_t = typename traits<ResidualTpl<PS>>::Meta_t;
        using ResidualModel_t = typename traits<ResidualMeta_t>::Model_t;
        using ResidualData_t = typename traits<ResidualMeta_t>::Data_t;

        using DimNH_t = typename traits<ResidualMeta_t>::DimNR_t;

        using H_t = Eigen::GMatrix<typename PS::VarScalar, DimNH_t::Value, 1, PS::Options>;
        using Hx_t = Eigen::GMatrix<typename PS::VarScalar, DimNH_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Hu_t = Eigen::GMatrix<typename PS::VarScalar, DimNH_t::Value, PS::DimNU_t::Value, PS::Options>;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl>
    struct traits<ConstraintDataResidualTpl<PhaseSpec, ResidualTpl>>
    {
        using SpecOfBaseClass = PhaseSpec;
        using Meta_t = ConstraintResidualTpl<PhaseSpec, ResidualTpl>;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl>
    struct traits<ConstraintModelResidualTpl<PhaseSpec, ResidualTpl>>
    {
        using SpecOfBaseClass = PhaseSpec;
        using Meta_t = ConstraintResidualTpl<PhaseSpec, ResidualTpl>;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl>
    struct ConstraintDataResidualTpl
        : public ConstraintDataBase<ConstraintDataResidualTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ConstraintResidualTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ConstraintDataBase<ConstraintDataResidualTpl<PS, ResidualTpl>, PS>;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(Meta_t);

        using ResidualMeta_t = typename traits<Meta_t>::ResidualMeta_t;
        using ResidualModel_t = typename traits<Meta_t>::ResidualModel_t;
        using ResidualData_t = typename traits<Meta_t>::ResidualData_t;

        DEFAULT_ACCESSOR(H_t, H);
        DEFAULT_ACCESSOR(Hx_t, Hx);
        DEFAULT_ACCESSOR(Hu_t, Hu);

        template <typename DataCollector>
        ConstraintDataResidualTpl(const Model_t &model, DataCollector *const collector)
            : residual(model.get_residual().createData(collector)),
              H(model.get_nh()),
              Hx(model.get_nh(), model.get_ps().get_ndx()),
              Hu(model.get_nh(), model.get_ps().get_nu())
        {
            H.setZero();
            Hx.setZero();
            Hu.setZero();
        }

        ResidualData_t residual;
        H_t H;
        Hx_t Hx;
        Hu_t Hu;

    }; // struct ConstraintDataResidualTpl

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl>
    class ConstraintModelResidualTpl
        : public ConstraintModelBase<ConstraintModelResidualTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ConstraintResidualTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ConstraintModelBase<ConstraintModelResidualTpl<PS, ResidualTpl>, PS>;

        using DimNH_t = typename traits<Meta_t>::DimNH_t;

        using ResidualMeta_t = typename traits<Meta_t>::ResidualMeta_t;
        using ResidualModel_t = typename traits<Meta_t>::ResidualModel_t;
        using ResidualData_t = typename traits<Meta_t>::ResidualData_t;

        ConstraintModelResidualTpl(const PS &ps, const ResidualModel_t &residual)
            : Base(ps, DimNH_t(residual.get_nr())),
              residual_(residual)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            residual_.calc(data.residual, x, u);
            data.H = data.residual.R;
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            residual_.calc(data.residual, x);
            data.H = data.residual.R;
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            residual_.calcDiff(data.residual, x, u);
            data.Hx = data.residual.Rx;
            data.Hu = data.residual.Ru;
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            residual_.calcDiff(data.residual, x);
            data.Hx = data.residual.Rx;
            data.Hu = data.residual.Ru;
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

        using Base::get_ps;

        using Base::get_nh;
        using Base::get_nh_dim;

    protected:
        ResidualModel_t residual_;

    }; // class ConstraintModelResidualTpl

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_residual_hpp__
