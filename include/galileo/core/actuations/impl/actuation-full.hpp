#ifndef __galileo_core_actuations_actuation_full_hpp__
#define __galileo_core_actuations_actuation_full_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

namespace galileo
{

    template <typename SystemSpec>
    class ActuationFullTpl;

    template <typename SystemSpec>
    struct traits<ActuationFullTpl<SystemSpec>>
    {
        using SS = SystemSpec;

        using Meta_t = ActuationFullTpl<SS>;
        using Model_t = ActuationModelFullTpl<SS>;
        using Data_t = ActuationDataFullTpl<SS>;

        using VectorNv_t = typename SS::VectorNv_t;
        using VectorNua_t = typename SS::VectorNua_t;
        using MatrixNvNdx_t = typename SS::MatrixNvNdx_t;
        using MatrixNvNua_t = typename SS::MatrixNvNua_t;
        using MatrixNuaNv_t = typename SS::MatrixNuaNv_t;
        using BoolArrayNv_t = Eigen::Array<bool, SS::NV, 1>;
    };

    template <typename SystemSpec>
    struct traits<ActuationModelFullTpl<SystemSpec>>
    {
        using Meta_t = ActuationFullTpl<SystemSpec>;
    };

    template <typename SystemSpec>
    struct traits<ActuationDataFullTpl<SystemSpec>>
    {
        using Meta_t = ActuationFullTpl<SystemSpec>;
    };

    template <typename SystemSpec>
    class ActuationDataFullTpl : public ActuationDataBase<ActuationDataFullTpl<SystemSpec>, SystemSpec>
    {
    public:
        using SS = SystemSpec;

        using Meta_t = ActuationFullTpl<SS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActuationDataBase<ActuationDataFullTpl<SS>, SS>;

        GALILEO_ACTUATION_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(VectorNv_t, tau);
        DEFAULT_ACCESSOR(VectorNua_t, u);
        DEFAULT_ACCESSOR(MatrixNvNdx_t, dtau_dx);
        DEFAULT_ACCESSOR(MatrixNvNua_t, dtau_du);
        DEFAULT_ACCESSOR(MatrixNuaNv_t, Mtau);
        DEFAULT_ACCESSOR(BoolArrayNv_t, tau_set);

        ActuationDataFullTpl(const Model_t &model)
            : tau(model.get_ss().get_nv()),
              u(model.get_ss().get_nua()),
              dtau_dx(model.get_ss().get_nv(), model.get_ss().get_ndx()),
              dtau_du(model.get_ss().get_nv(), model.get_ss().get_nua()),
              Mtau(model.get_ss().get_nua(), model.get_ss().get_nv()),
              tau_set(model.get_ss().get_nv())
        {
            tau.setZero();
            u.setZero();
            dtau_dx.setZero();
            dtau_du.setIdentity();
            Mtau.setIdentity();
            tau_set.setOnes();
        }

        VectorNv_t tau;
        VectorNua_t u;
        MatrixNvNdx_t dtau_dx;
        MatrixNvNua_t dtau_du;
        MatrixNuaNv_t Mtau;
        BoolArrayNv_t tau_set;
    };

    template <typename SystemSpec>
    class ActuationModelFullTpl : public ActuationModelBase<ActuationModelFullTpl<SystemSpec>, SystemSpec>
    {
    public:
        using SS = SystemSpec;

        using Meta_t = ActuationFullTpl<SS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActuationModelBase<ActuationModelFullTpl<SS>, SS>;

        using State_t = typename SS::State_t;

        ActuationModelFullTpl(const SS &ss, const State_t &state) : Base(ss, state) {}

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.tau = u;
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
            data.u = tau;
        }

        template <typename StateVectorType, typename ControlVectorType>
        void torqueTransform(Data_t &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            // has constant values which are set in createData
        }

        Data_t createData() const { return Data_t(*this); }

        using Base::get_ss;
        using Base::get_state;

    }; // class ActuationModelFullTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_full_hpp__
