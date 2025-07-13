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

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

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

        Data_t createData() const
        {
            return this->derived().createData();
        }

        const PS &get_ps() const
        {
            return ps_;
        }

        /**
         * @brief Return the dimension of the control space
         */
        const int get_nr() const
        {
            return this->derived().get_nr_impl();
        }

        const int get_nr_impl() const
        {
            if constexpr (DimNR_t::IsFixed)
            {
                return DimNR_t::Value;
            }
            else
            {
                return NR_dim_.value();
            }
        }

        const DimNR_t &NRDim() const
        {
            return NR_dim_;
        }

    protected:
        inline ActivationModelBase(const PS &ps, const DimNR_t &NR_dim) : ps_(ps), NR_dim_(NR_dim)
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

        const PS &ps_;
        const DimNR_t &NR_dim_;

    }; // class ActivationModelBase

} // namespace galileo

#endif // __galileo_core_activations_activation_model_base_hpp__
