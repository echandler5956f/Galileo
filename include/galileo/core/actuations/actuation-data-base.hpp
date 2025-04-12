#ifndef __galileo_core_actuations_actuation_data_base_hpp__
#define __galileo_core_actuations_actuation_data_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/core/actuations/actuation-model-base.hpp"

#include <array>

namespace galileo
{

    template <typename BasicSpec>
    struct ActuationDataTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using BS = BasicSpec;

        typename BS::VectorNv_t tau;
        typename BS::VectorNua_t u;
        typename BS::MatrixNvNdx_t dtau_dx;
        typename BS::MatrixNu_t dtau_du;
        typename BS::MatrixNuaNv_t Mtau;
        std::array<bool, BS::NV> tau_set;

    }; // struct ActuationDataTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_data_base_hpp__
