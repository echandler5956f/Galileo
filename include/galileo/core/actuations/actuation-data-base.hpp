#ifndef __galileo_core_actuations_actuation_data_base_hpp__
#define __galileo_core_actuations_actuation_data_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/core/actuations/actuation-model-base.hpp"

#include <array>

namespace galileo
{
    namespace core
    {

        template <typename PhaseSpec>
        struct ActuationDataTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            typename PS::VectorNv_t tau;
            typename PS::VectorNu_t u;
            typename PS::MatrixNvNdx_t dTaudX;
            typename PS::MatrixNvNu_t dTaudU;
            typename PS::MatrixNuNv_t Mtau;
            std::array<bool, PS::NV> tau_set;

        }; // struct ActuationDataTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_actuation_data_base_hpp__
