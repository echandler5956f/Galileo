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
            return ps_.get();
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
        inline ActivationModelBase(const PS &ps, const DimNR_t &nr_dim)
            : ps_(ps), nr_dim_(nr_dim)
        {
        }

        inline ActivationModelBase(const ActivationModelBase &clone)
        {
            *this = clone;
        }

        inline ActivationModelBase &operator=(const ActivationModelBase &clone)
        {
            ps_ = clone.ps_;
            nr_dim_ = clone.nr_dim_;
            return *this;
        }

        std::reference_wrapper<const PS> ps_;
        DimNR_t nr_dim_;

    }; // class ActivationModelBase

} // namespace galileo

#endif // __galileo_core_activations_activation_model_base_hpp__
