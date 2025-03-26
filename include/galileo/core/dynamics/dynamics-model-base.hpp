#ifndef __galileo_core_dynamics_dynamics_model_base_hpp__
#define __galileo_core_dynamics_dynamics_model_base_hpp__

#include "galileo/core/dynamics/dynamics-base.hpp"

#define GALILEO_DYNAMICS_BASIC_TYPEDEF(Dynamics)                                  \
    using VarScalar = typename traits<Dynamics>::VarScalar;                       \
    using NumScalar = typename traits<Dynamics>::NumScalar;                       \
    static constexpr int Options = traits<Dynamics>::Options;                     \
    using DynamicsModelDerived = typename traits<Dynamics>::DynamicsModelDerived; \
    using DynamicsDataDerived = typename traits<Dynamics>::DynamicsDataDerived;

#define GALILEO_DYNAMICS_CONSTANTS(Dynamics)

#define GALILEO_DYNAMICS_MODEL_TYPEDEF(Dynamics)        \
    using State_t = typename traits<Dynamics>::State_t; \
    using ActuationModel_t = typename traits<Dynamics>::ActuationModel_t;

#define GALILEO_DYNAMICS_DATA_TYPEDEF(Dynamics)                         \
    using ActuationData_t = typename traits<Dynamics>::ActuationData_t; \
    using F_t = typename traits<Dynamics>::F_t;                         \
    using Fx_t = typename traits<Dynamics>::Fx_t;                       \
    using Fu_t = typename traits<Dynamics>::Fu_t;

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

            template <typename StateVectorType, typename ControlVectorType>
            void calc(DynamicsDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(DynamicsDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
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

            State_t *state_;
            ActuationModel_t *actuation_;

        }; // class DynamicsModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_dynamics_dynamics_model_base_hpp__
