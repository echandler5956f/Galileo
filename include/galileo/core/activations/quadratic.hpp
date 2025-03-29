#ifndef __galileo_core_activations_quadratic_hpp__
#define __galileo_core_activations_quadratic_hpp__

#include "galileo/core/activations/activation-base.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationQuadraticTpl;

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ResidualModelTpl>
        struct traits<ActivationQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>>
        {
            using ResidualDerived = traits<_ResidualModelTpl<_VarScalar, _NumScalar, _Options>>::ResidualDerived;

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NR = traits<ResidualDerived>::NR;

            using ActivationDataDerived = ActivationDataQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>;
            using ActivationModelDerived = ActivationModelQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>;

            using A_t = VarScalar;
            using Ar_t = Eigen::Matrix<VarScalar, NR, 1, Options>;
            using Arr_t = Eigen::DiagonalMatrix<VarScalar, NR, Options>;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ResidualModelTpl>
        struct traits<ActivationDataQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>>
        {
            using ActivationDerived = ActivationQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>;
            using VarScalar = traits<ActivationDerived>::VarScalar;
            using NumScalar = traits<ActivationDerived>::NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ResidualModelTpl>
        struct traits<ActivationModelQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>>
        {
            using ActivationDerived = ActivationQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>;
            using VarScalar = traits<ActivationDerived>::VarScalar;
            using NumScalar = traits<ActivationDerived>::NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ResidualModelTpl>
        struct ActivationDataQuadraticTpl : public ActivationDataBase<ActivationDataQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActivationDerived = ActivationQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>;
            GALILEO_ACTIVATION_BASIC_TYPEDEF(ActivationDerived);
            GALILEO_ACTIVATION_CONSTANTS(ActivationDerived);
            GALILEO_ACTIVATION_DATA_TYPEDEF(ActivationDerived);

            A_t A;
            Ar_t Ar;
            Arr_t Arr;

        }; // class ActivationDataQuadraticTpl

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ResidualModelTpl>
        class ActivationModelQuadraticTpl : public ActivationModelBase<ActivationModelQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActivationDerived = ActivationQuadraticTpl<_VarScalar, _NumScalar, _Options, _ResidualModelTpl>;
            GALILEO_ACTIVATION_BASIC_TYPEDEF(ActivationDerived);
            GALILEO_ACTIVATION_CONSTANTS(ActivationDerived);
            GALILEO_ACTIVATION_MODEL_TYPEDEF(ActivationDerived);

            template <typename ResidualVectorType>
            void calc(ActivationDataDerived &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                data.A = VarScalar(0.5) * r.dot(r);
            }

            template <typename ResidualVectorType>
            void calcDiff(ActivationDataDerived &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                data.Ar = r;
                // The Hessian has constant values which were set in createData.
            }

        }; // class ActivationModelQuadraticTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_activations_quadratic_hpp__