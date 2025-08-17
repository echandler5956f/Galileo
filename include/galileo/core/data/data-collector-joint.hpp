#ifndef __galileo_core_data_data_collector_joint_hpp__
#define __galileo_core_data_data_collector_joint_hpp__

namespace galileo
{

    template <typename PhaseSpec>
    struct JointDataTpl
    {
        using PS = PhaseSpec;

        typename PS::VectorNua_t tau;        // Joint torques
        typename PS::VectorNv_t a;           // Joint accelerations
        typename PS::MatrixNuaNdx_t dtau_dx; // Torque derivatives w.r.t. state
        typename PS::MatrixNuaNu_t dtau_du;  // Torque derivatives w.r.t. control
        typename PS::MatrixNvNdx_t da_dx;    // Acceleration derivatives w.r.t. state
        typename PS::MatrixNvNu_t da_du;     // Acceleration derivatives w.r.t. control

        JointDataTpl(const PS &ps)
            : tau(ps.get_nua()),
              a(ps.get_nv()),
              dtau_dx(ps.get_nua(), ps.get_ndx()),
              dtau_du(ps.get_nua(), ps.get_nu()),
              da_dx(ps.get_nv(), ps.get_ndx()),
              da_du(ps.get_nv(), ps.get_nu())
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
        using PS = PhaseSpec;

        JointDataMixinTpl(std::shared_ptr<JointDataTpl<PS>> data) : joint(data) {}

        std::shared_ptr<JointDataTpl<PS>> joint;

    }; // struct JointDataMixinTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_joint_hpp__
