#ifndef __galileo_core_integrators_erk_base_hpp__
#define __galileo_core_integrators_erk_base_hpp__

#include "galileo/core/integrators/integrator-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        class IntegratorModelERKBase : IntegratorModelBase<IntegratorModelERKBase<Derived>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using IntegratorDerived = typename traits<Derived>::IntegratorDerived;
            GALILEO_INTEGRATOR_BASIC_TYPEDEF(IntegratorDerived);
            GALILEO_INTEGRATOR_CONSTANTS(IntegratorDerived);
            GALILEO_INTEGRATOR_MODEL_TYPEDEF(IntegratorDerived);

            template <typename StateVectorType, typename ControlVectorType>
            void calc(IntegratorDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x0,
                      const std::vector<Eigen::MatrixBase<ControlVectorType>> &us) const
                requires(IntegratorType == EXPLICIT_RK)
            {
                derived().calc(data, x0.derived(), us);
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(IntegratorDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x0,
                          const std::vector<Eigen::MatrixBase<ControlVectorType>> &us) const
                requires(IntegratorType == EXPLICIT_RK)
            {
                derived().calcDiff(data, x0.derived(), us);
            }

        protected:
            inline IntegratorModelERKBase()
            {
            }

            inline IntegratorModelERKBase(const IntegratorModelERKBase &clone)
            {
                *this = clone;
            }

            inline IntegratorModelERKBase &operator=(const IntegratorModelERKBase &clone)
            {
                return *this;
            }

        }; // class ERKBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_integrators_erk_base_hpp__
