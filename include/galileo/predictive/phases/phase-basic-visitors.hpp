#ifndef __galileo_predictive_phases_phase_basic_visitors_hpp__
#define __galileo_predictive_phases_phase_basic_visitors_hpp__

#include "galileo/predictive/phases/fwd.hpp"

#include "galileo/core/states/state-base.hpp"

namespace galileo
{

    namespace predictive
    {

        // These are basically all the visitors that we need in order to homogenously iterate over a set of heterogeneous phases

        // Phase model visitors

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename ControlVectorType>
        inline void segment_calc_zeroth_order(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &segment_index,
            const Eigen::MatrixBase<StateVectorType> &xs,
            const Eigen::MatrixBase<ControlVectorType> &us);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename ControlVectorType>
        inline void segment_calc_first_order(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &segment_index,
            const Eigen::MatrixBase<StateVectorType> &xs,
            const Eigen::MatrixBase<ControlVectorType> &us);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename ControlVectorType>
        inline void segment_quasi_static(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &segment_index,
            const Eigen::MatrixBase<StateVectorType> &xs,
            Eigen::MatrixBase<ControlVectorType> &us,
            const std::size_t &maxiter,
            const NumScalar &tol);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline NumScalar phase_period(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::VectorNX_t state_zero(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::VectorNX_t state_rand(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType1,
                  typename StateVectorType2,
                  typename StateTangentVectorType>
        inline void state_diff(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType1> &xs,
            const Eigen::MatrixBase<StateVectorType2> &xs_next,
            Eigen::MatrixBase<StateTangentVectorType> &dxout);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType,
                  typename StateVectorType2>
        inline void state_integrate(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx,
            Eigen::MatrixBase<StateVectorType2> &xout);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType1,
                  typename StateVectorType2,
                  typename JMatrix1,
                  typename JMatrix2>
        inline void state_jdiff(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType1> &x0,
            const Eigen::MatrixBase<StateVectorType2> &x1,
            Eigen::MatrixBase<JMatrix1> &Jfirst,
            Eigen::MatrixBase<JMatrix2> &Jsecond,
            const galileo::core::Jcomponent firstsecond = galileo::core::Jcomponent::both);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType,
                  typename JMatrix1,
                  typename JMatrix2>
        inline void state_jintegrate(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx,
            Eigen::MatrixBase<JMatrix1> &Jfirst,
            Eigen::MatrixBase<JMatrix2> &Jsecond,
            const galileo::core::Jcomponent firstsecond = galileo::core::Jcomponent::both,
            const galileo::core::AssignmentOp op = galileo::core::AssignmentOp::setto);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType,
                  typename JMatrix>
        inline void state_jintegrate_transport(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx,
            Eigen::MatrixBase<JMatrix> &Jin,
            const galileo::core::Jcomponent firstsecond);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType>
        inline typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::VectorNDX_t state_diff_dx(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x0,
            const Eigen::MatrixBase<StateVectorType> &x1);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType>
        inline typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::VectorNX_t state_integrate_x(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType1,
                  typename StateVectorType2>
        inline std::vector<typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::MatrixNDX_t> state_jdiff_Js(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType1> &x0,
            const Eigen::MatrixBase<StateVectorType2> &x1,
            const galileo::core::Jcomponent firstsecond = galileo::core::Jcomponent::both);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType>
        inline std::vector<typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::MatrixNDX_t> state_jintegrate_Js(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx,
            const galileo::core::Jcomponent firstsecond = galileo::core::Jcomponent::both);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int nx(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int nu(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int ndx(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int nh(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int ng(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int nc(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model);

        // Phase data visitors



    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_basic_visitors_hpp__