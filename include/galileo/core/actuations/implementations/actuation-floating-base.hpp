#ifndef __galileo_core_actuations_actuation_floating_base_hpp__
#define __galileo_core_actuations_actuation_floating_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

namespace galileo
{

    template <typename Spec>
    class ActuationFloatingBaseTpl;

    template <typename Spec>
    struct traits<ActuationFloatingBaseTpl<Spec>>
    {
        using SS = Spec;

        using Meta_t = ActuationFloatingBaseTpl<SS>;
        using Model_t = ActuationModelFloatingBaseTpl<SS>;
        using Data_t = ActuationDataFloatingBaseTpl<SS>;

        using VectorNv_t = typename SS::VectorNv_t;
        using VectorNua_t = typename SS::VectorNua_t;
        using MatrixNvNdx_t = typename SS::MatrixNvNdx_t;
        using MatrixNvNua_t = typename SS::MatrixNvNua_t;
        using MatrixNuaNv_t = typename SS::MatrixNuaNv_t;
        using BoolArrayNv_t = Eigen::Array<bool, SS::NV, 1>;
    };

    template <typename Spec>
    struct traits<ActuationModelFloatingBaseTpl<Spec>>
    {
        using Meta_t = ActuationFloatingBaseTpl<Spec>;
    };

    template <typename Spec>
    struct traits<ActuationDataFloatingBaseTpl<Spec>>
    {
        using Meta_t = ActuationFloatingBaseTpl<Spec>;
    };

    template <typename Spec>
    class ActuationDataFloatingBaseTpl : public ActuationDataBase<ActuationDataFloatingBaseTpl<Spec>, Spec>
    {
    public:
        using SS = Spec;

        using Meta_t = ActuationFloatingBaseTpl<SS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActuationDataBase<ActuationDataFloatingBaseTpl<SS>, SS>;

        GALILEO_ACTUATION_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(VectorNv_t, tau);
        DEFAULT_ACCESSOR(VectorNua_t, u);
        DEFAULT_ACCESSOR(MatrixNvNdx_t, dtau_dx);
        DEFAULT_ACCESSOR(MatrixNvNua_t, dtau_du);
        DEFAULT_ACCESSOR(MatrixNuaNv_t, Mtau);
        DEFAULT_ACCESSOR(BoolArrayNv_t, tau_set);

        ActuationDataFloatingBaseTpl(const Model_t &model)
            : tau(model.get_spec().get_nv()),
              u(model.get_spec().get_nua()),
              dtau_dx(model.get_spec().get_nv(), model.get_spec().get_ndx()),
              dtau_du(model.get_spec().get_nv(), model.get_spec().get_nua()),
              Mtau(model.get_spec().get_nua(), model.get_spec().get_nv()),
              tau_set(model.get_spec().get_nv())
        {
            tau.setZero();
            u.setZero();
            dtau_dx.setZero();
            dtau_du.setZero();
            Mtau.setZero();
            tau_set.setOnes();

            dtau_du.diagonal(-model.get_spec().get_nvb()).setOnes();
            Mtau.diagonal(model.get_spec().get_nvb()).setOnes();
            for (int i = 0; i < model.get_spec().get_nvb(); ++i) tau_set(i) = false;
        }

        VectorNv_t tau;
        VectorNua_t u;
        MatrixNvNdx_t dtau_dx;
        MatrixNvNua_t dtau_du;
        MatrixNuaNv_t Mtau;
        BoolArrayNv_t tau_set;
    };

    template <typename Spec>
    class ActuationModelFloatingBaseTpl : public ActuationModelBase<ActuationModelFloatingBaseTpl<Spec>, Spec>
    {
    public:
        using SS = Spec;

        using Meta_t = ActuationFloatingBaseTpl<SS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActuationModelBase<ActuationModelFloatingBaseTpl<SS>, SS>;

        using State_t = typename SS::State_t;

        ActuationModelFloatingBaseTpl(const SS &spec) : Base(spec) {}

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            tail(data.tau, get_spec().get_nua_dim()) = u;
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
            data.u = tail(tau, get_spec().get_nua_dim());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void torqueTransform(Data_t &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            // has constant values which are set in createData
        }

        Data_t createData() const { return Data_t(*this); }

        using Base::get_spec;

    }; // class ActuationModelFloatingBaseTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_floating_base_hpp__
