#ifndef __galileo_core_actuations_actuation_floating_base_hpp__
#define __galileo_core_actuations_actuation_floating_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/multibody/robot-spec.hpp"

namespace galileo
{

    template <typename RobotSpec>
    class ActuationFloatingBaseTpl;

    template <typename RobotSpec>
    struct traits<ActuationFloatingBaseTpl<RobotSpec>>
    {
        using RS = RobotSpec;

        using Meta_t = ActuationFloatingBaseTpl<RS>;
        using Model_t = ActuationModelFloatingBaseTpl<RS>;
        using Data_t = ActuationDataTpl<RS>;
    };

    template <typename RobotSpec>
    struct traits<ActuationModelFloatingBaseTpl<RobotSpec>>
    {
        using Meta_t = ActuationFloatingBaseTpl<RobotSpec>;
    };

    template <typename RobotSpec>
    class ActuationModelFloatingBaseTpl
        : public ActuationModelBase<ActuationModelFloatingBaseTpl<RobotSpec>, RobotSpec>
    {
    public:
        using RS = RobotSpec;

        using Meta_t = typename RS::ActuationMeta_t;
        using Model_t = typename RS::ActuationModel_t;
        using Data_t = typename RS::ActuationData_t;

        using Base = ActuationModelBase<ActuationModelFloatingBaseTpl<RS>, RS>;

        using State_t = typename RS::State_t;

        ActuationModelFloatingBaseTpl(const State_t &state)
            : Base(state)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            tail(data.tau, get_state().get_nua_dim()) = u;
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            // has constant values which are set in createData
        }

        template <typename StateVectorType, typename TauVectorType>
        void commands(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<TauVectorType> &tau) const
        {
            data.u = tail(tau, get_state().get_nua_dim());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void torqueTransform(Data_t &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            // has constant values which are set in createData
        }

        Data_t createData() const
        {
            Data_t data(*this);
            data.dtau_du.diagonal(-get_state().get_nvb()).setOnes();
            data.Mtau.diagonal(get_state().get_nvb()).setOnes();
            for (int i = 0; i < get_state().get_nvb(); ++i)
                data.tau_set(i) = false;
            return data;
        }

        using Base::get_state;

    }; // class ActuationModelFloatingBaseTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_floating_base_hpp__
