#ifndef __galileo_core_constraints_constraint_model_hpp__
#define __galileo_core_constraints_constraint_model_hpp__

#include "galileo/core/constraints/constraint-base.hpp"

namespace galileo
{
    namespace core
    {

        template <
            typename VarScalar,
            typename NumScalar,
            int Options,
            template <typename V, typename N, int O> class ResidualModelTpl,
            ConstraintType EqualityInequality>
        struct ConstraintResidualTpl;

        template <
            typename _VarScalar,
            typename _NumScalar,
            int _Options,
            template <typename V, typename N, int O> class _ResidualModelTpl,
            ConstraintType _EqualityInequality>
        struct traits<ConstraintResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>>
        {
            using ResidualTpl = traits<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>>::ResidualTpl;

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NX = traits<ResidualTpl>::NX;
            static constexpr int NU = traits<ResidualTpl>::NU;
            static constexpr int NH = constexpr(EqualityInequality == ConstraintType::Equality) ? traits<ResidualTpl>::NH : 0;
            static constexpr int NG = constexpr(EqualityInequality == ConstraintType::Inequality) ? traits<ResidualTpl>::NG : 0;

            using ConstraintDataDerived = ConstraintDataResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>;
            using ConstraintModelDerived = ConstraintModelResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>;

            using ResidualModel_t = traits<ResidualTpl>::ResidualModel_t;
            using ResidualData_t = traits<ResidualTpl>::ResidualData_t;

            using H_t = Eigen::Matrix<VarScalar, NH, 1, Options>;
            using Hx_t = Eigen::Matrix<VarScalar, NH, NX, Options>;
            using Hu_t = Eigen::Matrix<VarScalar, NH, NU, Options>;
            using G_t = Eigen::Matrix<VarScalar, NG, 1, Options>;
            using Gx_t = Eigen::Matrix<VarScalar, NG, NX, Options>;
            using Gu_t = Eigen::Matrix<VarScalar, NG, NU, Options>;
        };

        template <
            typename _VarScalar,
            typename _NumScalar,
            int _Options,
            template <typename V, typename N, int O> class _ResidualModelTpl,
            ConstraintType _EqualityInequality>
        struct traits<ConstraintDataResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>>
        {
            using ConstraintDerived = ConstraintResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>;
            using VarScalar = traits<ConstraintDerived>::VarScalar;
            using NumScalar = traits<ConstraintDerived>::NumScalar;
        };

        template <
            typename _VarScalar,
            typename _NumScalar,
            int _Options,
            template <typename V, typename N, int O> class _ResidualModelTpl,
            ConstraintType _EqualityInequality>
        struct traits<ConstraintModelResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>>
        {
            using ConstraintDerived = ConstraintResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>;
            using VarScalar = traits<ConstraintDerived>::VarScalar;
            using NumScalar = traits<ConstraintDerived>::NumScalar;
        };

        template <
            typename _VarScalar,
            typename _NumScalar,
            int _Options,
            template <typename V, typename N, int O> class _ResidualModelTpl,
            ConstraintType _EqualityInequality>
        struct ConstraintDataResidualTpl : ConstraintDataBase<ConstraintDataResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ConstraintDerived = ConstraintResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>;
            GALILEO_CONSTRAINT_BASIC_TYPEDEF(ConstraintDerived);
            GALILEO_CONSTRAINT_CONSTANTS(ConstraintDerived);
            GALILEO_CONSTRAINT_DATA_TYPEDEF(ConstraintDerived);

            ResidualData_t residual;
            H_t H;
            Hx_t Hx;
            Hu_t Hu;
            G_t G;
            Gx_t Gx;
            Gu_t Gu;

        }; // struct ConstraintDataResidualTpl

        template <
            typename _VarScalar,
            typename _NumScalar,
            int _Options,
            template <typename V, typename N, int O> class _ResidualModelTpl,
            ConstraintType _EqualityInequality>
        class ConstraintModelResidualTpl : public ConstraintModelBase<ConstraintModelResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ConstraintDerived = ConstraintResidualTpl<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>, _EqualityInequality>;
            GALILEO_CONSTRAINT_BASIC_TYPEDEF(ConstraintDerived);
            GALILEO_CONSTRAINT_CONSTANTS(ConstraintDerived);
            GALILEO_CONSTRAINT_MODEL_TYPEDEF(ConstraintDerived);

            ConstraintModelResidualTpl(const ResidualModel_t &residual)
                : residual_(residual)
            {
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calc(ConstraintDataDerived &data,
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
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                residual_.calcDiff(data.residual, x.derived(), u.derived());

                if constexpr (EqualityInequality == ConstraintType::Equality)
                    updateEqualityCalcDiff(data);
                else
                    updateInequalityCalcDiff(data);
            }

        protected:
            void updateEqualityCalc(ConstraintDataDerived &data) const
            {
                data.H = data.residual.r;
            }

            void updateInequalityCalc(ConstraintDataDerived &data) const
            {
                data.G = data.residual.r;
            }

            void updateEqualityCalcDiff(ConstraintDataDerived &data) const
            {
                data.Hx = data.residual.rx;
                data.Hu = data.residual.ru;
            }

            void updateInequalityCalcDiff(ConstraintDataDerived &data) const
            {
                data.Gx = data.residual.rx;
                data.Gu = data.residual.ru;
            }

            ResidualModel_t residual_;

        }; // class ConstraintModelResidualTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_model_hpp__
