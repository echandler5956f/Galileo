#ifndef __galileo_core_actuations_actuation_model_base_hpp__
#define __galileo_core_actuations_actuation_model_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

#define GALILEO_ACTUATIONS_BASIC_TYPEDEF(Actuations)                                    \
    using Scalar = typename traits<Actuations>::Scalar;                                 \
    using VarScalar = typename traits<Actuations>::VarScalar;                           \
    static constexpr int Options = traits<Actuations>::Options;                                \
    using ActuationsModelDerived = typename traits<Actuations>::ActuationsModelDerived; \
    using ActuationsDataDerived = typename traits<Actuations>::ActuationsDataDerived;

#define GALILEO_ACTUATIONS_CONSTANTS(Actuations) \

#define GALILEO_ACTUATIONS_MODEL_TYPEDEF(Actuations)

#define GALILEO_ACTUATIONS_DATA_TYPEDEF(Actuations) \

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        class ActuationsModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActuationsDerived = typename traits<Derived>::ActuationsDerived;
            GALILEO_ACTUATIONS_BASIC_TYPEDEF(ActuationsDerived);
            GALILEO_ACTUATIONS_CONSTANTS(ActuationsDerived);
            GALILEO_ACTUATIONS_MODEL_TYPEDEF(ActuationsDerived);

            template <typename ResidualVectorType>
            void calc(ActuationsDataDerived &data,
                      const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                derived().calc(data, r.derived());
            }

            template <typename ResidualVectorType>
            void calcDiff(ActuationsDataDerived &data,
                          const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                derived().calcDiff(data, r.derived());
            }

        protected:
            inline ActuationsModelBase()
            {
            }

            inline ActuationsModelBase(const ActuationsModelBase &clone)
            {
                *this = clone;
            }

            inline ActuationsModelBase &operator=(const ActuationsModelBase &clone)
            {
                return *this;
            }

        }; // class ActuationsModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_actuation_model_base_hpp__
