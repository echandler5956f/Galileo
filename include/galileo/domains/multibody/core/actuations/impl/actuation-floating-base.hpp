#ifndef __galileo_multibody_core_actuations_actuation_floating_base_hpp__
#define __galileo_multibody_core_actuations_actuation_floating_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/domains/multibody/core/actuations/fwd.hpp"

namespace galileo
{

    template <typename MultibodySpec>
    class ActuationFloatingBaseTpl;

    template <typename MultibodySpec>
    struct traits<ActuationFloatingBaseTpl<MultibodySpec>>
    {
        using SS = MultibodySpec;

        using Meta_t = ActuationFloatingBaseTpl<SS>;
        using Model_t = ActuationModelFloatingBaseTpl<SS>;
        using Data_t = ActuationDataFloatingBaseTpl<SS>;

        using VectorNv_t = ArenaMatrixTpl<typename SS::VectorNv_t>;
        using VectorNua_t = ArenaMatrixTpl<typename SS::VectorNua_t>;
        using MatrixNvNdx_t = ArenaMatrixTpl<typename SS::MatrixNvNdx_t>;
        using MatrixNvNua_t = ArenaMatrixTpl<typename SS::MatrixNvNua_t>;
        using MatrixNuaNv_t = ArenaMatrixTpl<typename SS::MatrixNuaNv_t>;
        using BoolArrayNv_t = ArenaMatrixTpl<Eigen::Array<bool, SS::NV, 1>>;
    };

    template <typename MultibodySpec>
    struct traits<ActuationModelFloatingBaseTpl<MultibodySpec>>
    {
        using Meta_t = ActuationFloatingBaseTpl<MultibodySpec>;
    };

    template <typename MultibodySpec>
    struct traits<ActuationDataFloatingBaseTpl<MultibodySpec>>
    {
        using Meta_t = ActuationFloatingBaseTpl<MultibodySpec>;
    };

    template <typename MultibodySpec>
    class ActuationDataFloatingBaseTpl
        : public ActuationDataBase<ActuationDataFloatingBaseTpl<MultibodySpec>, MultibodySpec>
    {
    public:
        using SS = MultibodySpec;

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

        ActuationDataFloatingBaseTpl(const Model_t &model, MemoryArena &arena)
            : tau(arena, model.get_ss().get_nv(), 1),
              u(arena, model.get_ss().get_nua(), 1),
              dtau_dx(arena, model.get_ss().get_nv(), model.get_ss().get_ndx()),
              dtau_du(arena, model.get_ss().get_nv(), model.get_ss().get_nua()),
              Mtau(arena, model.get_ss().get_nua(), model.get_ss().get_nv()),
              tau_set(arena, model.get_ss().get_nv(), 1)
        {
            tau.setZero();
            u.setZero();
            dtau_dx.setZero();
            dtau_du.setZero();
            Mtau.setZero();
            tau_set.setOnes();

            dtau_du.diagonal(-model.get_ss().get_nvb()).setOnes();
            Mtau.diagonal(model.get_ss().get_nvb()).setOnes();
            for (int i = 0; i < model.get_ss().get_nvb(); ++i) tau_set(i) = false;
        }

        VectorNv_t tau;
        VectorNua_t u;
        MatrixNvNdx_t dtau_dx;
        MatrixNvNua_t dtau_du;
        MatrixNuaNv_t Mtau;
        BoolArrayNv_t tau_set;

    }; // class ActuationDataFloatingBaseTpl

    template <typename MultibodySpec>
    class ActuationModelFloatingBaseTpl
        : public ActuationModelBase<ActuationModelFloatingBaseTpl<MultibodySpec>, MultibodySpec>
    {
    public:
        using SS = MultibodySpec;

        using Meta_t = ActuationFloatingBaseTpl<SS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActuationModelBase<ActuationModelFloatingBaseTpl<SS>, SS>;

        using State_t = typename SS::State_t;

        ActuationModelFloatingBaseTpl(const SS &ss, const State_t &state) : Base(ss, state) {}

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            tail(data.tau, get_ss().get_nua_dim()) = u;
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
            data.u = tail(tau, get_ss().get_nua_dim());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void torqueTransform(Data_t &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            // has constant values which are set in createData
        }

        Data_t createData(MemoryArena &arena) const { return Data_t(*this, arena); }

        using Base::get_ss;

    }; // class ActuationModelFloatingBaseTpl

} // namespace galileo

#endif // __galileo_multibody_core_actuations_actuation_floating_base_hpp__
