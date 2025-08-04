#ifndef __galileo_core_activations_activation_model_base_hpp__
#define __galileo_core_activations_activation_model_base_hpp__

#include "galileo/core/activations/activation-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ActivationModelBase
        : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        template <typename ResidualVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            this->derived().calc(data, r);
        }

        template <typename ResidualVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            this->derived().calcDiff(data, r);
        }

        Data_t createData() const
        {
            return this->derived().createData();
        }

        const DimNR_t &get_nr_dim() const
        {
            return nr_dim_;
        }

        int get_nr() const
        {
            return nr_dim_.value();
        }

    protected:
        inline ActivationModelBase(const DimNR_t &nr_dim)
            : nr_dim_(nr_dim)
        {
        }

        inline ActivationModelBase(const ActivationModelBase &clone)
            : nr_dim_(clone.nr_dim_)
        {
        }

        inline ActivationModelBase &operator=(const ActivationModelBase &clone)
        {
            nr_dim_ = clone.nr_dim_;
            return *this;
        }

        DimNR_t nr_dim_;

    }; // class ActivationModelBase

} // namespace galileo

#endif // __galileo_core_activations_activation_model_base_hpp__
