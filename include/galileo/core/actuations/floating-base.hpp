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
        using Data_t = ActuationDataTpl<BS>;
        using Model_t = ActuationModelFloatingBaseTpl<BS>;
    };

    template <typename BasicSpec>
    struct traits<ActuationModelFloatingBaseTpl<BasicSpec>>
    {
        using BS = BasicSpec;

        using Meta_t = ActuationFloatingBaseTpl<BS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
    };

    template <typename BasicSpec>
    class ActuationModelFloatingBaseTpl : public ActuationModelBase<ActuationModelFloatingBaseTpl<BasicSpec>, BasicSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using BS = BasicSpec;

        using ActuationDerived = ActuationFloatingBaseTpl<BS>;
        using ActuationDataDerived = typename traits<ActuationDerived>::ActuationDataDerived;
        using ActuationModelDerived = typename traits<ActuationDerived>::ActuationModelDerived;

        ActuationModelFloatingBaseTpl() {}

        template <typename StateVectorType, typename ControlVectorType>
        void calc(ActuationDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.tau.tail(BS::NUa) = u;
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(ActuationDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            // has constant values which are set in createData
        }

        template <typename StateVectorType, typename TauVectorType>
        void commands(ActuationDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<TauVectorType> &tau) const
        {
            data.u = tau.tail(BS::NUa);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void torqueTransform(ActuationDataDerived &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            // has constant values which are set in createData
        }

        ActuationDataDerived createData() const
        {
            ActuationDataDerived data = ActuationDataDerived();
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