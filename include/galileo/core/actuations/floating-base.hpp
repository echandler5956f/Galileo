#ifndef __galileo_core_actuations_floating - base_hpp__
#define __galileo_core_actuations_floating -base_hpp__

#include "galileo/core/actuations/actuation-model-base.hpp"

namespace galileo
{

    namespace core
    {

        template <typename BasicSpec>
        class ActuationModelFloatingBaseTpl : public ActuationModelBase<ActuationModelFloatingBaseTpl<BasicSpec>, BasicSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using BS = BasicSpec;

            template <typename StateVectorType, typename ControlVectorType>
            void calc(typename BS::ActuationData_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                data->tau.tail(BS::NUa) = u;
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(typename BS::ActuationData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                // has constant values which are set in createData
            }

            template <typename StateVectorType, typename TauVectorType>
            void commands(typename BS::ActuationData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<TauVectorType> &tau) const
            {
                data->u = tau.tail(BS::NUa);
            }

            template <typename StateVectorType, typename ControlVectorType>
            void torqueTransform(typename BS::ActuationData_t &data,
                                 const Eigen::MatrixBase<StateVectorType> &x,
                                 const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                // has constant values which are set in createData
            }

        }; // class ActuationModelFloatingBaseTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_floating-base_hpp__