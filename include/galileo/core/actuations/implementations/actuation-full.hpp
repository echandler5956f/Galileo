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
        using Data_t = ActuationDataTpl<RS>;
    };

    template <typename RobotSpec>
    struct traits<ActuationModelFullTpl<RobotSpec>>
    {
        using Meta_t = ActuationFullTpl<RobotSpec>;
    };

    template <typename RobotSpec>
    class ActuationModelFullTpl
        : public ActuationModelBase<ActuationModelFullTpl<RobotSpec>, RobotSpec>
    {
    public:
        using RS = RobotSpec;

        using Meta_t = typename RS::ActuationMeta_t;
        using Model_t = typename RS::ActuationModel_t;
        using Data_t = typename RS::ActuationData_t;

        using Base = ActuationModelBase<ActuationModelFullTpl<RS>, RS>;

        using State_t = typename RS::State_t;

        ActuationModelFullTpl(const State_t &state)
            : Base(state)
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
            Data_t data(*this);
            data.dtau_du.setIdentity();
            data.Mtau.setIdentity();
            return data;
        }

        using Base::get_state;

    }; // class ActuationModelFullTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_full_hpp__
