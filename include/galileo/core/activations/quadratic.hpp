#ifndef __galileo_core_activations_quadratic_hpp__
#define __galileo_core_activations_quadratic_hpp__

#include "galileo/core/activations/activation-base.hpp"

namespace galileo
{

    namespace core
    {

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationQuadraticTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct traits<ActivationQuadraticTpl<PhaseSpec, ResidualTpl>>
        {
            using PS = PhaseSpec;
            using ResidualMeta = traits<ResidualTpl<PS>>;
            using ResidualModel_t = typename traits<ResidualMeta>::ResidualModel_t;
            using ResidualData_t = typename traits<ResidualMeta>::ResidualData_t;

            static constexpr int NR = traits<ResidualMeta>::NR;

            using ActivationDataDerived = ActivationDataQuadraticTpl<PS, ResidualTpl>;
            using ActivationModelDerived = ActivationModelQuadraticTpl<PS, ResidualTpl>;

            using A_t = PS::VarScalar;
            using Ar_t = Eigen::Matrix<typename PS::VarScalar, NR, 1, PS::Options>;
            using Arr_t = Eigen::Matrix<typename PS::VarScalar, NR, NR, PS::Options>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct traits<ActivationDataQuadraticTpl<PhaseSpec, ResidualTpl>>
        {
            using PS = PhaseSpec;
            using ActivationDerived = ActivationQuadraticTpl<PS, ResidualTpl>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct traits<ActivationModelQuadraticTpl<PhaseSpec, ResidualTpl>>
        {
            using PS = PhaseSpec;
            using ActivationDerived = ActivationQuadraticTpl<PS, ResidualTpl>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataQuadraticTpl : public ActivationDataBase<ActivationDataQuadraticTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ActivationDerived = ActivationQuadraticTpl<PS, ResidualTpl>;
            GALILEO_ACTIVATION_DATA_TYPEDEF(ActivationDerived);

            DEFAULT_ACCESSOR(A_t, A);
            DEFAULT_ACCESSOR(Ar_t, Ar);
            DEFAULT_ACCESSOR(Arr_t, Arr);

            A_t A;
            Ar_t Ar;
            Arr_t Arr;

        }; // class ActivationDataQuadraticTpl

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelQuadraticTpl : public ActivationModelBase<ActivationModelQuadraticTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
            using ActivationDerived = ActivationQuadraticTpl<PS, ResidualTpl>;
            using ActivationDataDerived = typename traits<ActivationDerived>::ActivationDataDerived;
            using ActivationModelDerived = typename traits<ActivationDerived>::ActivationModelDerived;

            template <typename ResidualVectorType>
            void calc(ActivationDataDerived &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                data.A() = typename PS::VarScalar(0.5) * r.dot(r);
            }

            template <typename ResidualVectorType>
            void calcDiff(ActivationDataDerived &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                data.Ar() = r;
                // The Hessian has constant values which were set in createData.
            }

        }; // class ActivationModelQuadraticTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_activations_quadratic_hpp__