#ifndef __galileo_multibody_actuations_floating_base_hpp__
#define __galileo_multibody_actuations_floating_base_hpp__

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
        using RS = RobotSpec;

        using Meta_t = ActuationFloatingBaseTpl<RS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename RobotSpec>
    class ActuationModelFloatingBaseTpl : public ActuationModelBase<ActuationModelFloatingBaseTpl<RobotSpec>, RobotSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

        GALILEO_ROBOT_SPEC_EIGEN_TYPES_TYPEDEF(RS);

        using Meta_t = typename RS::ActuationMeta_t;
        using Model_t = typename RS::ActuationModel_t;
        using Data_t = typename RS::ActuationData_t;

        ActuationModelFloatingBaseTpl(std::shared_ptr<typename RS::State_t> state)
            : ActuationModelBase<ActuationModelFloatingBaseTpl<RS>, RS>(state) {}

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            tail(data.tau, this->NUaDim()) = u;
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
            data.u = tail(tau, this->NUaDim());
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
            data.dtau_du.diagonal(-this->get_state()->get_nvb()).setOnes();
            data.Mtau.diagonal(this->get_state()->get_nvb()).setOnes();
            for (int i = 0; i < this->get_state()->get_nvb(); ++i)
            {
                data.tau_set(i) = false;
            }
            return data;
        }

        using Base = ActuationModelBase<ActuationModelFloatingBaseTpl<RS>, RS>;

        using Base::get_state;

        using Base::get_nua;
        using Base::NUaDim;

    }; // class ActuationModelFloatingBaseTpl

} // namespace galileo

#endif // __galileo_multibody_actuations_floating_base_hpp__
