#ifndef __galileo_core_data_joint_hpp__
#define __galileo_core_data_joint_hpp__

#include "galileo/core/data/fwd.hpp"
#include "galileo/core/data/data-collector-base.hpp"
#include "galileo/core/data/data-collector-actuation.hpp"

namespace galileo
{

    namespace core
    {

        template <typename PhaseSpec>
        struct JointDataBaseTpl
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            typename PS::VectorNua_t tau;        // Joint torques
            typename PS::VectorNv_t a;           // Joint accelerations
            typename PS::MatrixNuaNdx_t dtau_dx; // Torque derivatives w.r.t. state
            typename PS::MatrixNu_t dtau_du;     // Torque derivatives w.r.t. control
            typename PS::MatrixNvNdx_t da_dx;    // Acceleration derivatives w.r.t. state
            typename PS::MatrixNvNu_t da_du;     // Acceleration derivatives w.r.t. control

        }; // struct JointDataBaseTpl

        // Joint data mixin
        template <typename Derived, typename PhaseSpec>
        struct JointDataMixinTpl
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            std::shared_ptr<JointDataBaseTpl<PhaseSpec>> joint;

            JointDataMixinTpl(std::shared_ptr<JointDataBaseTpl<PhaseSpec>> data)
                : joint(data) {}

        }; // struct JointDataMixinTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_data_joint_hpp__