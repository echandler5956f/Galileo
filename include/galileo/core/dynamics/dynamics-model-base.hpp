#ifndef __galileo_core_dynamics_dynamics_model_base_hpp__
#define __galileo_core_dynamics_dynamics_model_base_hpp__

#include "galileo/core/dynamics/dynamics-base.hpp"

#define GALILEO_DYNAMICS_BASIC_TYPEDEF(Dynamics)                                    \
    using Scalar = typename traits<Dynamics>::Scalar;                                 \
    using VarScalar = typename traits<Dynamics>::VarScalar;                           \
    static constexpr int Options = traits<Dynamics>::Options;                                \
    using DynamicsModelDerived = typename traits<Dynamics>::DynamicsModelDerived; \
    using DynamicsDataDerived = typename traits<Dynamics>::DynamicsDataDerived;

#define GALILEO_DYNAMICS_CONSTANTS(Dynamics) \

#define GALILEO_DYNAMICS_MODEL_TYPEDEF(Dynamics)

#define GALILEO_DYNAMICS_DATA_TYPEDEF(Dynamics) \

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        class DynamicsModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using DynamicsDerived = typename traits<Derived>::DynamicsDerived;
            GALILEO_DYNAMICS_BASIC_TYPEDEF(DynamicsDerived);
            GALILEO_DYNAMICS_CONSTANTS(DynamicsDerived);
            GALILEO_DYNAMICS_MODEL_TYPEDEF(DynamicsDerived);

            template <typename ResidualVectorType>
            void calc(DynamicsDataDerived &data,
                      const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                derived().calc(data, r.derived());
            }

            template <typename ResidualVectorType>
            void calcDiff(DynamicsDataDerived &data,
                          const Eigen::MatrixBase<ResidualVectorType> &r) const
            {
                derived().calcDiff(data, r.derived());
            }

        protected:
            inline DynamicsModelBase()
            {
            }

            inline DynamicsModelBase(const DynamicsModelBase &clone)
            {
                *this = clone;
            }

            inline DynamicsModelBase &operator=(const DynamicsModelBase &clone)
            {
                return *this;
            }

        }; // class DynamicsModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_dynamics_dynamics_model_base_hpp__
