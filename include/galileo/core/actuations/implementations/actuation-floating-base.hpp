#ifndef __galileo_core_actuations_actuation_floating_base_hpp__
#define __galileo_core_actuations_actuation_floating_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

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
        using Data_t = ActuationDataFloatingBaseTpl<RS>;

        using VectorNv_t = typename RS::VectorNv_t;
        using VectorNua_t = typename RS::VectorNua_t;
        using MatrixNvNdx_t = typename RS::MatrixNvNdx_t;
        using MatrixNvNua_t = typename RS::MatrixNvNua_t;
        using MatrixNuaNv_t = typename RS::MatrixNuaNv_t;
        using BoolArrayNv_t = Eigen::Array<bool, RS::NV, 1>;
    };

    template <typename RobotSpec>
    struct traits<ActuationModelFloatingBaseTpl<RobotSpec>>
    {
        using Meta_t = ActuationFloatingBaseTpl<RobotSpec>;
    };

    template <typename RobotSpec>
    struct traits<ActuationDataFloatingBaseTpl<RobotSpec>>
    {
        using Meta_t = ActuationFloatingBaseTpl<RobotSpec>;
    };

    template <typename RobotSpec>
    class ActuationDataFloatingBaseTpl
        : public ActuationDataBase<ActuationDataFloatingBaseTpl<RobotSpec>, RobotSpec>
    {
    public:
        using RS = RobotSpec;

        using Meta_t = ActuationFloatingBaseTpl<RS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActuationDataBase<ActuationDataFloatingBaseTpl<RS>, RS>;

        GALILEO_ACTUATION_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(VectorNv_t, tau);
        DEFAULT_ACCESSOR(VectorNua_t, u);
        DEFAULT_ACCESSOR(MatrixNvNdx_t, dtau_dx);
        DEFAULT_ACCESSOR(MatrixNvNua_t, dtau_du);
        DEFAULT_ACCESSOR(MatrixNuaNv_t, Mtau);
        DEFAULT_ACCESSOR(BoolArrayNv_t, tau_set);

        ActuationDataFloatingBaseTpl(const Model_t &model)
            : tau(model.get_state().get_nv()),
              u(model.get_state().get_nua()),
              dtau_dx(model.get_state().get_nv(),
                      model.get_state().get_ndx()),
              dtau_du(model.get_state().get_nv(),
                      model.get_state().get_nua()),
              Mtau(model.get_state().get_nua(),
                   model.get_state().get_nv()),
              tau_set(model.get_state().get_nv())
        {
            tau.setZero();
            u.setZero();
            dtau_dx.setZero();
            dtau_du.setZero();
            Mtau.setZero();
            tau_set.setOnes();

            dtau_du.diagonal(-model.get_state().get_nvb()).setOnes();
            Mtau.diagonal(model.get_state().get_nvb()).setOnes();
            for (int i = 0; i < model.get_state().get_nvb(); ++i)
                tau_set(i) = false;
        }

        VectorNv_t tau;
        VectorNua_t u;
        MatrixNvNdx_t dtau_dx;
        MatrixNvNua_t dtau_du;
        MatrixNuaNv_t Mtau;
        BoolArrayNv_t tau_set;
    };

    template <typename RobotSpec>
    class ActuationModelFloatingBaseTpl
        : public ActuationModelBase<ActuationModelFloatingBaseTpl<RobotSpec>, RobotSpec>
    {
    public:
        using RS = RobotSpec;

        using Meta_t = ActuationFloatingBaseTpl<RS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
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
            return Data_t(*this);
        }

        using Base::get_state;

    }; // class ActuationModelFloatingBaseTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_floating_base_hpp__
