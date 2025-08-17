#ifndef __galileo_multibody_impulses_impulse_visitors_hpp__
#define __galileo_multibody_impulses_impulse_visitors_hpp__

#include "galileo/multibody/impulses/fwd.hpp"

namespace galileo
{

    // These are basically all the visitors that we need in order to homogenously iterate over a set of heterogeneous impulses

    // Impulse model visitors

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename StateVectorType>
    inline void impulse_calc_zeroth_order(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                          ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                                          const Eigen::MatrixBase<StateVectorType> &x);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename StateVectorType>
    inline void impulse_calc_first_order(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                         ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                                         const Eigen::MatrixBase<StateVectorType> &x);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename ForceVectorType>
    inline void impulse_update_force(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                     ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                                     const Eigen::MatrixBase<ForceVectorType> &force);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename MatrixNcNdxType>
    inline void impulse_update_force_diff(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                          ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                                          const Eigen::MatrixBase<MatrixNcNdxType> &df_dx);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline void impulse_set_zero_force(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                       ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline void impulse_set_zero_force_diff(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                            ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> impulse_create_data(
        const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
        typename PhaseSpec::RobotData_t *const robot);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t impulse_get_id(
        const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline void impulse_set_id(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                               const typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t &id);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t impulse_get_type(
        const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline void impulse_set_type(
        const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
        const typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t &type);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline int impulse_get_nc(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model);

    // Impulse data visitors

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::RobotDataPointer_t impulse_robot_data(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t impulse_frame(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t impulse_type_data(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::SE3_t impulse_jMf(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNv_t impulse_Jc(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::Force_t impulse_f(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::Force_t impulse_fext(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNdx_t impulse_df_dx(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNu_t impulse_df_du(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::ActionMatrix_t impulse_fXj(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNv_t impulse_dv0_dq(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNv_t impulse_dtau_dq(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data);

} // namespace galileo

#endif // __galileo_multibody_impulses_impulse_visitors_hpp__
