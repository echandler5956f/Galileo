#ifndef __galileo_core_activations_activation_model_base_hpp__
#define __galileo_core_activations_activation_model_base_hpp__

#include "galileo/core/activations/activation-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ActivationModelBase : internal::CRTP<ActivationModelBase<Derived, PhaseSpec>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using ActivationDerived = typename traits<Derived>::ActivationDerived;
        using ActivationDataDerived = typename traits<ActivationDerived>::ActivationDataDerived;
        using ActivationModelDerived = typename traits<ActivationDerived>::ActivationModelDerived;

        template <typename ResidualVectorType>
        void calc(ActivationDataDerived &data,
                  const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            this->derived().calc(data, r.derived());
        }

        template <typename ResidualVectorType>
        void calcDiff(ActivationDataDerived &data,
                      const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            this->derived().calcDiff(data, r.derived());
        }

        template <typename DataCollector>
        ActivationDataDerived createData(DataCollector *const collector)
        {
            return this->derived().createData(collector);
        }

        int nr() const
        {
            return this->derived().nr_impl();
        }

        int nr_impl() const
        {
            return traits<ActivationDerived>::NR;
        }

    protected:
        inline ActivationModelBase()
        {
        }

        inline ActivationModelBase(const ActivationModelBase &clone)
        {
            *this = clone;
        }

        inline ActivationModelBase &operator=(const ActivationModelBase &clone)
        {
            return *this;
        }

    }; // class ActivationModelBase

} // namespace galileo

#endif // __galileo_core_activations_activation_model_base_hpp__
