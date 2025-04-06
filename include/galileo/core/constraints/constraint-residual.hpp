#ifndef __galileo_core_constraints_constraint_model_hpp__
#define __galileo_core_constraints_constraint_model_hpp__

#include "galileo/core/constraints/constraint-base.hpp"

namespace galileo
{
    namespace core
    {

        template <
            typename PhaseSpec,
            template <typename PS> class ResidualModelTpl,
            ConstraintType EqualityInequality>
        struct ConstraintResidualTpl;

        template <
            typename PhaseSpec,
            template <typename PS> class ResidualTpl,
            ConstraintType EqualityInequality_>
        struct traits<ConstraintResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality_>>
        {
            using PS = PhaseSpec;
            using ResidualMeta = traits<ResidualTpl<PS>>;
            using ResidualModel_t = traits<ResidualMeta>::ResidualModel_t;
            using ResidualData_t = traits<ResidualMeta>::ResidualData_t;

            using EqualityInequality = EqualityInequality_;
            static constexpr int NH = constexpr(EqualityInequality == typename ConstraintType::Equality) ? traits<ResidualMeta>::NR : 0;
            static constexpr int NG = constexpr(EqualityInequality == typename ConstraintType::Inequality) ? traits<ResidualMeta>::NR : 0;

            using ConstraintDataDerived = ConstraintDataResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality>;
            using ConstraintModelDerived = ConstraintModelResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality>;

            using H_t = Eigen::Matrix<typename PS::VarScalar, NH, 1, PS::Options>;
            using Hx_t = Eigen::Matrix<typename PS::VarScalar, NH, PS::NDX, PS::Options>;
            using Hu_t = Eigen::Matrix<typename PS::VarScalar, NH, PS::NU, PS::Options>;
            using G_t = Eigen::Matrix<typename PS::VarScalar, NG, 1, PS::Options>;
            using Gx_t = Eigen::Matrix<typename PS::VarScalar, NG, PS::NDX, PS::Options>;
            using Gu_t = Eigen::Matrix<typename PS::VarScalar, NG, PS::NU, PS::Options>;
        };

        template <
            typename PhaseSpec,
            template <typename PS> class ResidualTpl,
            ConstraintType EqualityInequality>
        struct traits<ConstraintDataResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality>>
        {
            using PS = PhaseSpec;
            using ConstraintDerived = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality>;
        };

        template <
            typename PhaseSpec,
            template <typename PS> class ResidualTpl,
            ConstraintType EqualityInequality>
        struct traits<ConstraintModelResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality>>
        {
            using PS = PhaseSpec;
            using ConstraintDerived = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality>;
        };

        template <
            typename PhaseSpec,
            template <typename PS> class ResidualTpl,
            ConstraintType EqualityInequality>
        struct ConstraintDataResidualTpl : ConstraintDataBase<ConstraintDataResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality>, PhaseSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ConstraintDerived = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality>;
            GALILEO_CONSTRAINT_DATA_TYPEDEF(ConstraintDerived);

            using ResidualData_t = traits<ResidualTpl<PS>>::ResidualData_t;

            DEFAULT_ACCESSOR(H_t, H);
            DEFAULT_ACCESSOR(Hx_t, Hx);
            DEFAULT_ACCESSOR(Hu_t, Hu);
            DEFAULT_ACCESSOR(G_t, G);
            DEFAULT_ACCESSOR(Gx_t, Gx);
            DEFAULT_ACCESSOR(Gu_t, Gu);

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
        class ConstraintModelResidualTpl : public ConstraintModelBase<ConstraintModelResidualTpl<PhaseSpec, ResidualTpl, EqualityInequality_>, PhaseSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
            static constexpr ConstraintType EqualityInequality = EqualityInequality_;
            using Constraint_t = ConstraintResidualTpl<PS, ResidualTpl, EqualityInequality>;
            using ConstraintModel_t = traits<Constraint_t>::ConstraintModelDerived;
            using ConstraintData_t = traits<Constraint_t>::ConstraintDataDerived;

            using BoundVector_t = Eigen::Matrix<typename PS::NumScalar, traits<Constraint_t>::NG, 1, PS::Options>;

            ConstraintModelResidualTpl(const ResidualModel_t &residual)
                : residual_(residual)
            {
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calc(ConstraintData_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                residual_.calc(data.residual, x.derived(), u.derived());

                if constexpr (EqualityInequality == ConstraintType::Equality)
                    updateEqualityCalc(data);
                else
                    updateInequalityCalc(data);
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(ConstraintData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                residual_.calcDiff(data.residual, x.derived(), u.derived());

                if constexpr (EqualityInequality == ConstraintType::Equality)
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

            const BoundVector_t &lb() const
            {
                return lb_;
            }

            const BoundVector_t &ub() const
            {
                return ub_;
            }

        protected:
            void updateEqualityCalc(ConstraintData_t &data) const
            {
                data.H = data.residual.r;
            }

            void updateInequalityCalc(ConstraintData_t &data) const
            {
                data.G = data.residual.r;
            }

            void updateEqualityCalcDiff(ConstraintData_t &data) const
            {
                data.Hx = data.residual.rx;
                data.Hu = data.residual.ru;
            }

            void updateInequalityCalcDiff(ConstraintData_t &data) const
            {
                data.Gx = data.residual.rx;
                data.Gu = data.residual.ru;
            }

            ResidualModel_t residual_;
            BoundVector_t lb_;
            BoundVector_t ub_;

        }; // class ConstraintModelResidualTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_model_hpp__
