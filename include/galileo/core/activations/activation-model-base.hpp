#ifndef __galileo_core_activations_activation_model_base_hpp__
#define __galileo_core_activations_activation_model_base_hpp__

#include "galileo/core/activations/activation-base.hpp"

#define GALILEO_ACTIVATION_BASIC_TYPEDEF(Activation)                                    \
    using Scalar = typename traits<Activation>::Scalar;                                 \
    using VarScalar = typename traits<Activation>::VarScalar;                           \
    static constexpr int Options = traits<Activation>::Options;                                \
    using ActivationModelDerived = typename traits<Activation>::ActivationModelDerived; \
    using ActivationDataDerived = typename traits<Activation>::ActivationDataDerived;

#define GALILEO_ACTIVATION_CONSTANTS(Activation) \
    static constexpr int NR = traits<Activation>::NR;

#define GALILEO_ACTIVATION_MODEL_TYPEDEF(Activation)

#define GALILEO_ACTIVATION_DATA_TYPEDEF(Activation) \
    using A_t = typename traits<Activation>::A_t;   \
    using Ar_t = typename traits<Activation>::Ar_t; \
    using Arr_t = typename traits<Activation>::Arr_t;

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        class ActivationModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActivationDerived = typename traits<Derived>::ActivationDerived;
            GALILEO_ACTIVATION_BASIC_TYPEDEF(ActivationDerived);
            GALILEO_ACTIVATION_CONSTANTS(ActivationDerived);
            GALILEO_ACTIVATION_MODEL_TYPEDEF(ActivationDerived);

            template <typename ResidualVectorType>
            void calc(ActivationDataDerived &data,
                      const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                derived().calc(data, r.derived());
            }

            template <typename ResidualVectorType>
            void calcDiff(ActivationDataDerived &data,
                          const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                derived().calcDiff(data, r.derived());
            }

        protected:
            inline ActivationModelBase()
            {
            }

            inline ActivationModelBase(const ActivationModelBase &clone)
            {
                *this = clone;
            }

            inline ActivationModelBase &operator=(const ActivationModelBase &clone)
            {
                return *this;
            }

        }; // class ActivationModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_activations_activation_model_base_hpp__
