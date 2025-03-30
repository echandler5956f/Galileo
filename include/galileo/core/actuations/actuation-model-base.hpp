#ifndef __galileo_core_actuations_actuation_model_base_hpp__
#define __galileo_core_actuations_actuation_model_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

#define GALILEO_ACTUATIONS_BASIC_TYPEDEF(Actuation)                                  \
    using Scalar = typename traits<Actuation>::Scalar;                               \
    using VarScalar = typename traits<Actuation>::VarScalar;                         \
    static constexpr int Options = traits<Actuation>::Options;                       \
    using ActuationModelDerived = typename traits<Actuation>::ActuationModelDerived; \
    using ActuationDataDerived = typename traits<Actuation>::ActuationDataDerived;

#define GALILEO_ACTUATIONS_CONSTANTS(Actuation)

#define GALILEO_ACTUATIONS_MODEL_TYPEDEF(Actuation)

#define GALILEO_ACTUATIONS_DATA_TYPEDEF(Actuation)                 \
    using VectorTau_t = typename traits<Actuation>::VectorTau_t;   \
    using VectorU_t = typename traits<Actuation>::VectorU_t;       \
    using MatrixTauX_t = typename traits<Actuation>::MatrixTauX_t; \
    using MatrixTauU_t = typename traits<Actuation>::MatrixTauU_t; \
    using MatrixMtau_t = typename traits<Actuation>::MatrixMtau_t;

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename PhaseSpec>
        class ActuationModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActuationDerived = typename traits<Derived>::ActuationDerived;
            GALILEO_ACTUATIONS_BASIC_TYPEDEF(ActuationDerived);
            GALILEO_ACTUATIONS_CONSTANTS(ActuationDerived);
            GALILEO_ACTUATIONS_MODEL_TYPEDEF(ActuationDerived);

            template <typename StateVectorType, typename ControlVectorType>
            void calc(ActuationDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(ActuationDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename TauVectorType>
            void commands(ActuationDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<TauVectorType> &tau) const
            {
                derived().commands(data, x.derived(), tau.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void torqueTransform(ActuationDataDerived &data,
                                 const Eigen::MatrixBase<StateVectorType> &x,
                                 const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().torqueTransform(data, x.derived(), u.derived());
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
