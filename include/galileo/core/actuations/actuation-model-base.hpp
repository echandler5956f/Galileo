#ifndef __galileo_core_actuations_actuation_model_base_hpp__
#define __galileo_core_actuations_actuation_model_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename BasicSpec>
        class ActuationModelBase : internal::CRTP<ActuationModelBase<Derived, BasicSpec>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using BS = BasicSpec;

            template <typename StateVectorType, typename ControlVectorType>
            void calc(typename BS::ActuationData_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                this->derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(typename BS::ActuationData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                this->derived().calcDiff(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename TauVectorType>
            void commands(typename BS::ActuationData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<TauVectorType> &tau) const
            {
                this->derived().commands(data, x.derived(), tau.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void torqueTransform(typename BS::ActuationData_t &data,
                                 const Eigen::MatrixBase<StateVectorType> &x,
                                 const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                this->derived().torqueTransform(data, x.derived(), u.derived());
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
