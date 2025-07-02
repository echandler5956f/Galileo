#ifndef __galileo_core_activations_activation_model_base_hpp__
#define __galileo_core_activations_activation_model_base_hpp__

#include "galileo/core/activations/activation-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ActivationModelBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        template <typename ResidualVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            this->derived().calc(data, r.derived());
        }

        template <typename ResidualVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            this->derived().calcDiff(data, r.derived());
        }

        Data_t createData()
        {
            return this->derived().createData();
        }

        int nr() const
        {
            return this->derived().nr_impl();
        }

        int nr_impl() const
        {
            return traits<Meta_t>::NR;
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
