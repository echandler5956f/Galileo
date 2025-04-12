#ifndef __galileo_core_actuations_floating_base_hpp__
#define __galileo_core_actuations_floating_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

namespace galileo
{

    template <typename BasicSpec>
    class ActuationFloatingBaseTpl;

    template <typename BasicSpec>
    struct traits<ActuationFloatingBaseTpl<BasicSpec>>
    {
        using BS = BasicSpec;

        using Meta_t = ActuationFloatingBaseTpl<BS>;
        using Model_t = ActuationModelFloatingBaseTpl<BS>;
        using Data_t = ActuationDataTpl<BS>;
    };

    template <typename BasicSpec>
    struct traits<ActuationModelFloatingBaseTpl<BasicSpec>>
    {
        using BS = BasicSpec;

        using Meta_t = ActuationFloatingBaseTpl<BS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename BasicSpec>
    class ActuationModelFloatingBaseTpl : public ActuationModelBase<ActuationModelFloatingBaseTpl<BasicSpec>, BasicSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using BS = BasicSpec;

        using Meta_t = typename BS::ActuationMeta_t;
        using Model_t = typename BS::ActuationModel_t;
        using Data_t = typename BS::ActuationData_t;

        ActuationModelFloatingBaseTpl() {}

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.tau.tail(BS::NUa) = u;
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
            data.u = tau.tail(BS::NUa);
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
            data.dtau_du.diagonal(-BS::NVb).setOnes();
            data.Mtau.diagonal(BS::NVb).setOnes();
            for (std::size_t i = 0; i < BS::NVb; ++i)
            {
                data.tau_set[i] = false;
            }
            return data;
        }

    }; // class ActuationModelFloatingBaseTpl

} // namespace galileo

#endif // __galileo_core_actuations_floating_base_hpp__