#ifndef __galileo_core_activations_activation_data_base_hpp__
#define __galileo_core_activations_activation_data_base_hpp__

#include "galileo/core/activations/activation-base.hpp"
#include "galileo/core/activations/activation-model-base.hpp"

// We use traits rather than PhaseSpec,
// because each activation model has its own NR
#define GALILEO_ACTIVATION_DATA_TYPEDEF(Activation)   \
    using A_t = typename traits<Activation>::A_t;     \
    using Ar_t = typename traits<Activation>::Ar_t;   \
    using Arr_t = typename traits<Activation>::Arr_t; \
    using Arr_diag_t = typename traits<Activation>::Arr_diag_t;

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct ActivationDataBase : internal::CRTP<ActivationDataBase<Derived, PhaseSpec>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using ActivationDerived = typename traits<Derived>::ActivationDerived;
        using ActivationDataDerived = typename traits<ActivationDerived>::ActivationDataDerived;
        using ActivationModelDerived = typename traits<ActivationDerived>::ActivationModelDerived;

        GALILEO_ACTIVATION_DATA_TYPEDEF(ActivationDerived);

        FORWARD_ACCESSOR(A_t, A);
        FORWARD_ACCESSOR(Ar_t, Ar);
        FORWARD_ACCESSOR(Arr_t, Arr);

        static Arr_diag_t getHessianMatrix(const ActivationDataDerived &data)
        {
            return data.Arr.diagonal().asDiagonal();
        }

        static void setHessianMatrix(ActivationDataDerived &data, const Arr_t &Arr)
        {
            data.Arr.diagonal() = Arr.diagonal();
        }

    protected:
        inline ActivationDataBase()
        {
        }

        inline ActivationDataBase(const ActivationDataBase &clone)
        {
            *this = clone;
        }

        inline ActivationDataBase &operator=(const ActivationDataBase &clone)
        {
            return *this;
        }

    }; // struct ActivationDataBase

} // namespace galileo

#endif // __galileo_core_activations_activation_data_base_hpp__
