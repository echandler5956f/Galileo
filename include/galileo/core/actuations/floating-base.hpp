#ifndef __galileo_core_actuations_floating - base_hpp__
#define __galileo_core_actuations_floating -base_hpp__

#include "galileo/core/actuations/actuation-model-base.hpp"

namespace galileo
{

    namespace core
    {

        template <typename PhaseSpec, int _NFbv>
        class ActuationModelFloatingBaseTpl : public ActuationModelBase<ActuationModelFloatingBaseTpl<PhaseSpec, _NFbv>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
            static constexpr int NFbv = _NFbv;
            static constexpr int Ntau = PS::NV - NFbv;

            template <typename StateVectorType, typename ControlVectorType>
            void calc(typename PS::ActuationData_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                data->tau.tail(Ntau) = u;
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(typename PS::ActuationData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                // has constant values which are set in createData
            }

            template <typename StateVectorType, typename TauVectorType>
            void commands(typename PS::ActuationData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<TauVectorType> &tau) const
            {
                data->u.tail(Ntau) = tau;
            }

            template <typename StateVectorType, typename ControlVectorType>
            void torqueTransform(typename PS::ActuationData_t &data,
                                 const Eigen::MatrixBase<StateVectorType> &x,
                                 const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                // has constant values which are set in createData
            }

        }; // class ActuationModelFloatingBaseTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_floating-base_hpp__