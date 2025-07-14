#ifndef __galileo_core_constraints_constraint_residual_hpp__
#define __galileo_core_constraints_constraint_residual_hpp__

#include "galileo/core/constraints/constraint-base.hpp"

namespace galileo
{

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        ConstraintType EqualityInequality>
    struct ConstraintResidualTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        ConstraintType EqualityInequality_>
    struct traits<ConstraintResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality_>>
    {
        using PS = PhaseSpec;

        using Meta_t = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality_>;
        using Model_t = ConstraintModelResidualTpl<PS, ResidualTpl, EqualityInequality_>;
        using Data_t = ConstraintDataResidualTpl<PS, ResidualTpl, EqualityInequality_>;

        using ResidualMeta_t = typename traits<ResidualTpl<PS>>::Meta_t;
        using ResidualModel_t = typename traits<ResidualMeta_t>::Model_t;
        using ResidualData_t = typename traits<ResidualMeta_t>::Data_t;

        static constexpr ConstraintType EqualityInequality = EqualityInequality_;
        using NHDim_t = is_equality_v<EqualityInequality> ? traits<ResidualMeta_t>::DimNR_t : DimensionTpl<0>;
        using NGDim_t = is_inequality_v<EqualityInequality> ? traits<ResidualMeta_t>::DimNR_t : DimensionTpl<0>;

        static constexpr int NH = NHDim_t::Value;
        static constexpr int NG = NGDim_t::Value;

        using H_t = Eigen::GMatrix<typename PS::VarScalar, NH, 1, PS::Options>;
        using Hx_t = Eigen::GMatrix<typename PS::VarScalar, NH, PS::DimNDX_t::Value, PS::Options>;
        using Hu_t = Eigen::GMatrix<typename PS::VarScalar, NH, PS::DimNU_t::Value, PS::Options>;
        using G_t = Eigen::GMatrix<typename PS::VarScalar, NG, 1, PS::Options>;
        using Gx_t = Eigen::GMatrix<typename PS::VarScalar, NG, PS::DimNDX_t::Value, PS::Options>;
        using Gu_t = Eigen::GMatrix<typename PS::VarScalar, NG, PS::DimNU_t::Value, PS::Options>;

        using BoundVector_t = Eigen::GMatrix<typename PS::NumScalar, NG, 1, PS::Options>;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        ConstraintType EqualityInequality>
    struct traits<ConstraintDataResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality>>
    {
        using PS = PhaseSpec;

        using Meta_t = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        ConstraintType EqualityInequality>
    struct traits<ConstraintModelResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality>>
    {
        using PS = PhaseSpec;

        using Meta_t = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        ConstraintType EqualityInequality>
    struct ConstraintDataResidualTpl
        : public ConstraintDataBase<ConstraintDataResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(Meta_t);

        using ResidualMeta_t = typename traits<Meta_t>::ResidualMeta_t;
        using ResidualModel_t = typename traits<Meta_t>::ResidualModel_t;
        using ResidualData_t = typename traits<Meta_t>::ResidualData_t;

        DEFAULT_ACCESSOR(H_t, H);
        DEFAULT_ACCESSOR(Hx_t, Hx);
        DEFAULT_ACCESSOR(Hu_t, Hu);
        DEFAULT_ACCESSOR(G_t, G);
        DEFAULT_ACCESSOR(Gx_t, Gx);
        DEFAULT_ACCESSOR(Gu_t, Gu);

        ConstraintDataResidualTpl(const Model_t &model)
            : residual(model.get_residual().createData(),
                       H(model.get_nh()),
                       Hx(model.get_nh(), model.get_ps().DimNDX_t::Value),
                       Hu(model.get_nh(), model.get_ps().DimNU_t::Value),
                       G(model.get_ng()),
                       Gx(model.get_ng(), model.get_ps().DimNDX_t::Value),
                       Gu(model.get_ng(), model.get_ps().DimNU_t::Value))
        {
            H.setZero();
            Hx.setZero();
            Hu.setZero();
            G.setZero();
            Gx.setZero();
            Gu.setZero();
        }

        ResidualData_t residual;
        H_t H;
        Hx_t Hx;
        Hu_t Hu;
        G_t G;
        Gx_t Gx;
        Gu_t Gu;

    }; // struct ConstraintDataResidualTpl

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        ConstraintType EqualityInequality_>
    class ConstraintModelResidualTpl
        : public ConstraintModelBase<ConstraintModelResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality_>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality_>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Base = ConstraintModelBase<ConstraintModelResidualTpl<PS, ResidualTpl, EqualityInequality_>, PS>;

        using DimNH_t = typename traits<Meta_t>::DimNH_t;
        using DimNG_t = typename traits<Meta_t>::DimNG_t;

        using ResidualMeta_t = typename traits<Meta_t>::ResidualMeta_t;
        using ResidualModel_t = typename traits<Meta_t>::ResidualModel_t;
        using ResidualData_t = typename traits<Meta_t>::ResidualData_t;

        static constexpr ConstraintType EqualityInequality = traits<Meta_t>::EqualityInequality;
        using BoundVector_t = typename traits<Meta_t>::BoundVector_t;

        ConstraintModelResidualTpl(const PS &ps, const ResidualModel_t &residual)
            : Base(ps, DimNH_t(is_equality_v<EqualityInequality> ? residual.get_nr() : 0),
                   DimNG_t(is_inequality_v<EqualityInequality> ? residual.get_nr() : 0)),
              residual_(residual)
        {
        }

        ConstraintModelResidualTpl(const PS &ps, const ResidualModel_t &residual,
                                   const BoundVector_t &lb, const BoundVector_t &ub)
            : Base(ps, DimNH_t(is_equality_v<EqualityInequality> ? residual.get_nr() : 0),
                   DimNG_t(is_inequality_v<EqualityInequality> ? residual.get_nr() : 0)),
              residual_(residual),
              lb_(lb),
              ub_(ub)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            residual_.calc(data.residual, x.derived(), u.derived());

            if constexpr (is_equality_v<EqualityInequality>)
                updateEqualityCalc(data);
            else
                updateInequalityCalc(data);
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            residual_.calc(data.residual, x.derived());

            if constexpr (is_equality_v<EqualityInequality>)
                updateEqualityCalc(data);
            else
                updateInequalityCalc(data);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            residual_.calcDiff(data.residual, x.derived(), u.derived());

            if constexpr (is_equality_v<EqualityInequality>)
                updateEqualityCalcDiff(data);
            else
                updateInequalityCalcDiff(data);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            residual_.calcDiff(data.residual, x.derived());

            if constexpr (is_equality_v<EqualityInequality>)
                updateEqualityCalcDiff(data);
            else
                updateInequalityCalcDiff(data);
        }

        template <typename LowerBoundType, typename UpperBoundType>
        void updateBounds(const Eigen::MatrixBase<LowerBoundType> &lb,
                          const Eigen::MatrixBase<UpperBoundType> &ub)
        {
            lb_ = lb.derived();
            ub_ = ub.derived();
        }

        const ResidualModel_t &get_residual() const
        {
            return residual_;
        }

        const BoundVector_t &get_lb() const
        {
            return lb_;
        }

        const BoundVector_t &get_ub() const
        {
            return ub_;
        }

        using Base::get_ps;

        using Base::get_nh;
        using Base::NHDim;

        using Base::get_ng;
        using Base::NGDim;

    protected:
        void updateEqualityCalc(Data_t &data) const
        {
            data.H = data.residual.R;
        }

        void updateInequalityCalc(Data_t &data) const
        {
            data.G = data.residual.R;
        }

        void updateEqualityCalcDiff(Data_t &data) const
        {
            data.Hx = data.residual.Rx;
            data.Hu = data.residual.Ru;
        }

        void updateInequalityCalcDiff(Data_t &data) const
        {
            data.Gx = data.residual.Rx;
            data.Gu = data.residual.Ru;
        }

        ResidualModel_t residual_;
        BoundVector_t lb_;
        BoundVector_t ub_;

    }; // class ConstraintModelResidualTpl

} // namespace galileo

#endif // __galileo_core_constraints_constraint_residual_hpp__
