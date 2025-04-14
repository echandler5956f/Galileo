#ifndef __galileo_core_actuations_actuation_data_base_hpp__
#define __galileo_core_actuations_actuation_data_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/core/actuations/actuation-model-base.hpp"

#include <array>

namespace galileo
{

    template <typename RobotSpec>
    struct ActuationDataTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

        using Meta_t = typename RS::ActuationMeta_t;
        using Model_t = typename RS::ActuationModel_t;
        using Data_t = typename RS::ActuationData_t;

        typename RS::VectorNv_t tau;
        typename RS::VectorNua_t u;
        typename RS::MatrixNvNdx_t dtau_dx;
        typename RS::MatrixNvNua_t dtau_du;
        typename RS::MatrixNuaNv_t Mtau;
        std::array<bool, RS::NV> tau_set;

    }; // struct ActuationDataTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_data_base_hpp__
