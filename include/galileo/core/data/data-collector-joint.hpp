#ifndef __galileo_core_data_data_collector_joint_hpp__
#define __galileo_core_data_data_collector_joint_hpp__

namespace galileo
{

    template <typename PhaseSpec>
    struct JointDataTpl
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        typename PS::VectorNua_t tau;        // Joint torques
        typename PS::VectorNv_t a;           // Joint accelerations
        typename PS::MatrixNuaNdx_t dtau_dx; // Torque derivatives w.r.t. state
        typename PS::MatrixNuaNu_t dtau_du;  // Torque derivatives w.r.t. control
        typename PS::MatrixNvNdx_t da_dx;    // Acceleration derivatives w.r.t. state
        typename PS::MatrixNvNu_t da_du;     // Acceleration derivatives w.r.t. control

        JointDataTpl(int nu)
            : tau(PS::NUa),
              a(PS::NV),
              dtau_dx(PS::NUa, PS::NDX),
              dtau_du(PS::NUa, nu),
              da_dx(PS::NV, PS::NDX),
              da_du(PS::NV, nu)
        {
            tau.setZero();
            a.setZero();
            dtau_dx.setZero();
            dtau_du.setZero();
            da_dx.setZero();
            da_du.setZero();
        }

    }; // struct JointDataTpl

    // Joint data mixin
    template <typename Derived, typename PhaseSpec>
    struct JointDataMixinTpl
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        JointDataTpl<PS> *joint;

        JointDataMixinTpl(JointDataTpl<PS> *data)
            : joint(data) {}

    }; // struct JointDataMixinTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_joint_hpp__
