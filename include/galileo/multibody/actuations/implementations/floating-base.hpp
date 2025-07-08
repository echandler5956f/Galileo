#ifndef __galileo_multibody_actuations_floating_base_hpp__
#define __galileo_multibody_actuations_floating_base_hpp__

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

        using Meta_t = typename RS::ActuationMeta_t;
        using Model_t = typename RS::ActuationModel_t;
        using Data_t = typename RS::ActuationData_t;

        ActuationModelFloatingBaseTpl() {}

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.tau.tail(RS::NUa) = u;
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
            data.u = tau.tail(RS::NUa);
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
            Data_t data = Data_t();
            data.dtau_du.diagonal(-RS::NVb).setOnes();
            data.Mtau.diagonal(RS::NVb).setOnes();
            for (std::size_t i = 0; i < RS::NVb; ++i)
            {
                data.tau_set[i] = false;
            }
            return data;
        }

    }; // class ActuationModelFloatingBaseTpl

} // namespace galileo

#endif // __galileo_multibody_actuations_floating_base_hpp__