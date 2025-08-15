#ifndef __galileo_core_actuations_actuation_full_hpp__
#define __galileo_core_actuations_actuation_full_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

namespace galileo
{

    template <typename RobotSpec>
    class ActuationFullTpl;

    template <typename RobotSpec>
    struct traits<ActuationFullTpl<RobotSpec>>
    {
        using RS = RobotSpec;

        using Meta_t = ActuationFullTpl<RS>;
        using Model_t = ActuationModelFullTpl<RS>;
        using Data_t = ActuationDataFullTpl<RS>;

        using VectorNv_t = typename RS::VectorNv_t;
        using VectorNua_t = typename RS::VectorNua_t;
        using MatrixNvNdx_t = typename RS::MatrixNvNdx_t;
        using MatrixNvNua_t = typename RS::MatrixNvNua_t;
        using MatrixNuaNv_t = typename RS::MatrixNuaNv_t;
        using BoolArrayNv_t = Eigen::Array<bool, RS::NV, 1>;
    };

    template <typename RobotSpec>
    struct traits<ActuationModelFullTpl<RobotSpec>>
    {
        using Meta_t = ActuationFullTpl<RobotSpec>;
    };

    template <typename RobotSpec>
    struct traits<ActuationDataFullTpl<RobotSpec>>
    {
        using Meta_t = ActuationFullTpl<RobotSpec>;
    };

     template <typename RobotSpec>
    class ActuationDataFullTpl
        : public ActuationDataBase<ActuationDataFullTpl<RobotSpec>, RobotSpec>
    {
    public:
        using RS = RobotSpec;

        using Meta_t = ActuationFullTpl<RS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActuationDataBase<ActuationDataFullTpl<RS>, RS>;

        GALILEO_ACTUATION_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(VectorNv_t, tau);
        DEFAULT_ACCESSOR(VectorNua_t, u);
        DEFAULT_ACCESSOR(MatrixNvNdx_t, dtau_dx);
        DEFAULT_ACCESSOR(MatrixNvNua_t, dtau_du);
        DEFAULT_ACCESSOR(MatrixNuaNv_t, Mtau);
        DEFAULT_ACCESSOR(BoolArrayNv_t, tau_set);

        ActuationDataFullTpl(const Model_t &model)
            : tau(model.get_rs().get_nv()),
              u(model.get_rs().get_nua()),
              dtau_dx(model.get_rs().get_nv(),
                      model.get_rs().get_ndx()),
              dtau_du(model.get_rs().get_nv(),
                      model.get_rs().get_nua()),
              Mtau(model.get_rs().get_nua(),
                   model.get_rs().get_nv()),
              tau_set(model.get_rs().get_nv())
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

    template <typename RobotSpec>
    class ActuationModelFullTpl
        : public ActuationModelBase<ActuationModelFullTpl<RobotSpec>, RobotSpec>
    {
    public:
        using RS = RobotSpec;

        using Meta_t = ActuationFullTpl<RS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActuationModelBase<ActuationModelFullTpl<RS>, RS>;

        using State_t = typename RS::State_t;

        ActuationModelFullTpl(const State_t &state)
            : Base(state.get_rs())
        {
        }

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

        Data_t createData() const
        {
            return Data_t(*this);
        }

        using Base::get_rs;

    }; // class ActuationModelFullTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_full_hpp__
