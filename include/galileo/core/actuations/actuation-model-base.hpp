#ifndef __galileo_core_actuations_actuation_model_base_hpp__
#define __galileo_core_actuations_actuation_model_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

#define GALILEO_ACTUATIONS_BASIC_TYPEDEF(Actuation)                                    \
    using Scalar = typename traits<Actuation>::Scalar;                                 \
    using VarScalar = typename traits<Actuation>::VarScalar;                           \
    static constexpr int Options = traits<Actuation>::Options;                         \
    using ActuationModelDerived = typename traits<Actuation>::ActuationModelDerived; \
    using ActuationDataDerived = typename traits<Actuation>::ActuationDataDerived;

#define GALILEO_ACTUATIONS_CONSTANTS(Actuation)

#define GALILEO_ACTUATIONS_MODEL_TYPEDEF(Actuation)

#define GALILEO_ACTUATIONS_DATA_TYPEDEF(Actuation)

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        class ActuationModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActuationDerived = typename traits<Derived>::ActuationDerived;
            GALILEO_ACTUATIONS_BASIC_TYPEDEF(ActuationDerived);
            GALILEO_ACTUATIONS_CONSTANTS(ActuationDerived);
            GALILEO_ACTUATIONS_MODEL_TYPEDEF(ActuationDerived);

            template <typename ResidualVectorType>
            void calc(ActuationDataDerived &data,
                      const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                derived().calc(data, r.derived());
            }

            template <typename ResidualVectorType>
            void calcDiff(ActuationDataDerived &data,
                          const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                derived().calcDiff(data, r.derived());
            }

        protected:
            inline ActuationModelBase()
            {
            }

            inline ActuationModelBase(const ActuationModelBase &clone)
            {
                *this = clone;
            }

            inline ActuationModelBase &operator=(const ActuationModelBase &clone)
            {
                return *this;
            }

        }; // class ActuationModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_actuation_model_base_hpp__
