#ifndef __galileo_core_actuations_actuation_data_base_hpp__
#define __galileo_core_actuations_actuation_data_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/core/actuations/actuation-model-base.hpp"
#include "galileo/multibody/robot-spec.hpp"

#include <array>

namespace galileo
{

    template <typename RobotSpec>
    struct ActuationDataTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(RS);

        using Meta_t = typename RS::ActuationMeta_t;
        using Model_t = typename RS::ActuationModel_t;
        using Data_t = typename RS::ActuationData_t;

        ActuationDataTpl(const Model_t &model) : tau(model.get_state()->get_nv()),
                                                 u(model.get_state()->get_nua()),
                                                 dtau_dx(model.get_state()->get_nv(),
                                                         model.get_state()->get_ndx()),
                                                 dtau_du(model.get_state()->get_nv(),
                                                         model.get_state()->get_nua()),
                                                 Mtau(model.get_state()->get_nua(),
                                                      model.get_state()->get_nv()),
                                                 tau_set(model.get_state()->get_nv())
        {
            tau.setZero();
            u.setZero();
            dtau_dx.setZero();
            dtau_du.setZero();
            Mtau.setZero();
            tau_set.setOnes();
        }

        VectorNv_t tau;
        VectorNua_t u;
        MatrixNvNdx_t dtau_dx;
        MatrixNvNua_t dtau_du;
        MatrixNuaNv_t Mtau;
        Eigen::Array<bool, RS::NV, 1> tau_set;

    }; // struct ActuationDataTpl

} // namespace galileo

#endif // __galileo_core_actuations_actuation_data_base_hpp__
