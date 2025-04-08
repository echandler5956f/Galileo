#ifndef __galileo_predictive_segments_segment_erk_234_hpp__
#define __galileo_predictive_segments_segment_erk_234_hpp__

#include "galileo/predictive/segments/segment-erk-base.hpp"

namespace galileo
{
    namespace predictive
    {

        enum class RKType
        {
            RK2 = 2,
            RK3 = 3,
            RK4 = 4
        };

        template <typename PhaseSpec, RKType _RKType>
        struct SegmentERKData234Tpl : public SegmentERKDataBase<SegmentERKData234Tpl<PhaseSpec, _RKType>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            typename PS::NodeDataVector_t nodes;
            typename PS::ControlParamDataVector_t controls;

            typename PS::VarScalarArray_t integral;
            typename PS::VectorNdx_t dx;
            typename PS::VectorNdxArray_t ki;
            typename PS::VectorNxArray_t y;
            typename PS::VectorNuArray_t us;
            typename PS::VectorNdxArray_t dx_rk;

            typename PS::MatrixNdxArray_t dki_dx;
            typename PS::MatrixNdxNwArray_t dki_dw;

            typename PS::MatrixNdxArray_t dyi_dx;
            typename PS::MatrixNdxNwArray_t dyi_dw;

            typename PS::VectorNdxArray_t dli_dx;
            typename PS::VectorNwArray_t dli_dw;

            typename PS::MatrixNdxArray_t ddli_ddx;
            typename PS::MatrixNuArray_t ddli_ddu;
            typename PS::MatrixNwArray_t ddli_ddw;
            typename PS::MatrixNdxNuArray_t ddli_dxdu;
            typename PS::MatrixNdxNwArray_t ddli_dxdw;
            typename PS::MatrixNuNwArray_t ddli_dudw;

            typename PS::MatrixNwArray_t Lww_partialx;
            typename PS::MatrixNdxNwArray_t Lxw_i;
            typename PS::MatrixNdxArray_t Lxx_partialx;
            typename PS::MatrixNdxNwArray_t Lxx_partialw;

            typename PS::XNext_t XNext;
            typename PS::Fx_t Fx;
            typename PS::Fw_t Fw;
            typename PS::L_t L;
            typename PS::Lx_t Lx;
            typename PS::Lw_t Lw;
            typename PS::Lxx_t Lxx;
            typename PS::Lxw_t Lxw;
            typename PS::Lww_t Lww;
            typename PS::H_t H;
            typename PS::Hx_t Hx;
            typename PS::Hw_t Hw;
            typename PS::G_t G;
            typename PS::Gx_t Gx;
            typename PS::Gw_t Gw;

        }; // struct SegmentERKData234Tpl

        template <typename PhaseSpec, RKType _RKType>
        class SegmentERKModel234Tpl : public SegmentERKModelBase<SegmentERKModel234Tpl<PhaseSpec, _RKType>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            template <typename StateVectorType, typename ControlParamVectorType>
            void calc(typename PS::SegmentData_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlParamVectorType> &ws) const
            {
                typename PS::NodeData_t &k0_data = data.nodes[0];
                typename PS::ControlParamData_t &w0_data = data.controls[0];

                control_.calc(w0_data, timings_[0], ws);
                data.us[0] = w0_data.u;
                node_.calc(k0_data, x, data.us[0]);
                data.y[0] = x;
                data.ki[0].head(PS::NV) = data.y[0].tail(PS::NV);
                data.ki[0].tail(PS::NV) = k0_data.XAcc;
                data.integral[0] = k0_data.L;
                for (std::size_t i = 1; i < PS::NStages; ++i)
                {
                    typename PS::NodeData_t &ki_data = data.nodes[i];
                    typename PS::ControlParamData_t &wi_data = data.controls[i];
                    data.dx_rk[i].noalias() = period_ * timings_[i] * data.ki[i - 1];
                    state_.integrate(x, data.dx_rk[i], data.y[i]);
                    control_.calc(wi_data, timings_[i], ws);
                    data.us[i] = wi_data.u;
                    node_.calc(ki_data, data.y[i], data.us[i]);
                    data.ki[i].head(PS::NV) = data.y[i].tail(PS::NV);
                    data.ki[i].tail(PS::NV) = ki_data.XAcc;
                    data.integral[i] = ki_data.L;
                }

                if (PS::NStages == 2)
                {
                    data.dx = data.ki[1] * period_;
                    data.L = data.integral[1] * period_;
                }
                else if (PS::NStages == 3)
                {
                    data.dx = (data.ki[0] + typename PS::NumScalar(3.) * data.ki[2]) * period_ /
                              typename PS::NumScalar(4.);
                    data.L = (data.integral[0] + typename PS::NumScalar(3.) * data.integral[2]) *
                             period_ / typename PS::NumScalar(4.);
                }
                else
                {
                    data.dx = (data.ki[0] + typename PS::NumScalar(2.) * data.ki[1] +
                               typename PS::NumScalar(2.) * data.ki[2] + data.ki[3]) *
                              period_ / typename PS::NumScalar(6.);
                    data.L = (data.integral[0] + typename PS::NumScalar(2.) * data.integral[1] +
                              typename PS::NumScalar(2.) * data.integral[2] + data.integral[3]) *
                             period_ / typename PS::NumScalar(6.);
                }
                state_.integrate(x, data.dx, data.XNext);
                data.H = k0_data.H;
                data.G = k0_data.G;
                if (with_cost_residual_)
                {
                    data.r = k0_data.r;
                }
            }

            template <typename StateVectorType, typename ControlParamVectorType>
            void calcDiff(typename PS::SegmentData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlParamVectorType> &ws) const
            {

                for (std::size_t i = 0; i < PS::NStages; ++i)
                {
                    node_.calcDiff(data.nodes[i], data.y[i], data.us[i]);
                }

                typename PS::NodeData_t &k0_data = data.nodes[0];
                typename PS::ControlParamData_t &w0_data = data.controls[0];
                data.dki_dx[0].bottomRows(PS::NV) = k0_data.Fx;
                control_.multiplyByJacobian(
                    w0_data, k0_data.Fu,
                    data.dki_dw[0].bottomRows(PS::NV)); // dki_dw = dki_du * du_dw

                data.dli_dx[0] = k0_data.Lx;
                control_.multiplyJacobianTransposeBy(
                    w0_data, k0_data.Lu,
                    data.dli_dw[0]); // dli_dw = dli_du * du_dw

                data.ddli_ddx[0] = k0_data.Lxx;
                data.ddli_ddu[0] = k0_data.Luu;
                control_.multiplyByJacobian(
                    w0_data, data.ddli_ddu[0],
                    data.ddli_dudw[0]); // ddli_dudw = ddli_ddu * du_dw
                control_.multiplyJacobianTransposeBy(
                    w0_data, data.ddli_dudw[0],
                    data.ddli_ddw[0]); // ddli_ddw = du_dw.T * ddli_dudw
                data.ddli_dxdu[0] = k0_data.Lxu;
                control_.multiplyByJacobian(
                    w0_data, data.ddli_dxdu[0],
                    data.ddli_dxdw[0]); // ddli_dxdw = ddli_dxdu * du_dw

                for (std::size_t i = 1; i < PS::NStages; ++i)
                {
                    typename PS::NodeData_t &ki_data = data.nodes[i];
                    typename PS::ControlParamData_t &wi_data = data.controls[i];
                    data.dyi_dx[i].noalias() = data.dki_dx[i - 1] * timings_[i] * period_;
                    data.dyi_dw[i].noalias() = data.dki_dw[i - 1] * timings_[i] * period_;
                    state_.JintegrateTransport(x, data.dx_rk[i], data.dyi_dx[i], second);
                    state_.Jintegrate(x, data.dx_rk[i], data.dyi_dx[i], data.dyi_dx[i], first, addto);
                    state_.JintegrateTransport(x, data.dx_rk[i], data.dyi_dw[i], second); // dyi_dw = Jintegrate * dyi_dw

                    // Sparse matrix-matrix multiplication for computing:
                    Eigen::Block<typename PS::MatrixNv_t> dkvi_dq = data.dki_dx[i].bottomLeftCorner(PS::NV, PS::NV);
                    Eigen::Block<typename PS::MatrixNv_t> dkvi_dv = data.dki_dx[i].bottomRightCorner(PS::NV, PS::NV);
                    Eigen::Block<typename PS::MatrixNvNw_t> dkqi_dw = data.dki_dw[i].topLeftCorner(PS::NV, PS::NW);
                    Eigen::Block<typename PS::MatrixNvNw_t> dkvi_dw = data.dki_dw[i].bottomLeftCorner(PS::NV, PS::NW);
                    const Eigen::Block<typename PS::MatrixNv_t> dki_dqi = ki_data.Fx.bottomLeftCorner(PS::NV, PS::NV);
                    const Eigen::Block<typename PS::MatrixNv_t> dki_dvi = ki_data.Fx.bottomRightCorner(PS::NV, PS::NV);

                    const Eigen::Block<typename PS::MatrixNv_t> dqi_dq = data.dyi_dx[i].topLeftCorner(PS::NV, PS::NV);
                    const Eigen::Block<typename PS::MatrixNv_t> dqi_dv = data.dyi_dx[i].topRightCorner(PS::NV, PS::NV);
                    const Eigen::Block<typename PS::MatrixNv_t> dvi_dq = data.dyi_dx[i].bottomLeftCorner(PS::NV, PS::NV);
                    const Eigen::Block<typename PS::MatrixNv_t> dvi_dv = data.dyi_dx[i].bottomRightCorner(PS::NV, PS::NV);
                    const Eigen::Block<typename PS::MatrixNvNW_t> dqi_dw = data.dyi_dw[i].topLeftCorner(PS::NV, PS::NW);
                    const Eigen::Block<typename PS::MatrixNvNW_t> dvi_dw = data.dyi_dw[i].bottomLeftCorner(PS::NV, PS::NW);
                    //   i. data.dki_dx[i].noalias() = data.dki_dy[i] * data.dyi_dx[i], where dki_dy
                    //   is ki_data.Fx
                    data.dki_dx[i].topRows(PS::NV) = data.dyi_dx[i].bottomRows(PS::NV);
                    dkvi_dq.noalias() = dki_dqi * dqi_dq;
                    if (i == 1)
                    {
                        dkvi_dv = period_ / typename PS::NumScalar(2.) * dki_dqi;
                    }
                    else
                    {
                        dkvi_dv.noalias() = dki_dqi * dqi_dv;
                    }
                    dkvi_dq.noalias() += dki_dvi * dvi_dq;
                    dkvi_dv.noalias() += dki_dvi * dvi_dv;
                    //  ii. data.dki_dw[i].noalias() = data.dki_dy[i] * data.dyi_dw[i], where dki_dy
                    //  is ki_data.Fx
                    dkqi_dw = dvi_dw;
                    dkvi_dw.noalias() = dki_dqi * dqi_dw;
                    dkvi_dw.noalias() += dki_dvi * dvi_dw;

                    control_.multiplyByJacobian(wi_data, ki_data.Fu,
                                                data.dki_dw[i].bottomRows(PS::NV),
                                                addto); // dfi_dw = dki_du * du_dw

                    data.dli_dx[i].noalias() = ki_data.Lx.transpose() * data.dyi_dx[i];
                    control_.multiplyJacobianTransposeBy(wi_data, ki_data.Lu,
                                                         data.dli_dw[i]); // dli_dw = Lu * du_dw
                    data.dli_dw[i].noalias() += ki_data.Lx.transpose() * data.dyi_dw[i];

                    data.Lxx_partialx[i].noalias() = ki_data.Lxx * data.dyi_dx[i];
                    data.ddli_ddx[i].noalias() = data.dyi_dx[i].transpose() * data.Lxx_partialx[i];

                    control_.multiplyByJacobian(wi_data, ki_data.Lxu,
                                                data.Lxw_i[i]); // Lxw = Lxu * du_dw
                    data.Lww_partialx[i].noalias() = data.Lxw_i[i].transpose() * data.dyi_dw[i];
                    data.Lxx_partialw[i].noalias() = ki_data.Lxx * data.dyi_dw[i];
                    control_.multiplyByJacobian(
                        wi_data, ki_data.Luu,
                        data.ddli_dudw[i]); // ddli_dudw = ddli_ddu * du_dw
                    control_.multiplyJacobianTransposeBy(
                        wi_data, data.ddli_dudw[i],
                        data.ddli_ddw[i]); // ddli_ddw = du_dw.T * ddli_dudw
                    data.ddli_ddw[i].noalias() += data.Lww_partialx[i].transpose() +
                                                  data.Lww_partialx[i] +
                                                  data.dyi_dw[i].transpose() * data.Lxx_partialw[i];

                    data.ddli_dxdu[i].noalias() = data.dyi_dx[i].transpose() * ki_data.Lxu;
                    control_.multiplyByJacobian(
                        wi_data, data.ddli_dxdu[i],
                        data.ddli_dxdw[i]); // ddli_dxdw = ddli_dxdu * du_dw
                    data.ddli_dxdw[i].noalias() += data.dyi_dx[i].transpose() * data.Lxx_partialw[i];
                }

                if (PS::NStages == 2)
                {
                    data.Fx.noalias() = period_ * data.dki_dx[1];
                    data.Fw.noalias() = period_ * data.dki_dw[1];
                    data.Lx.noalias() = period_ * data.dli_dx[1];
                    data.Lw.noalias() = period_ * data.dli_dw[1];
                    data.Lxx.noalias() = period_ * data.ddli_ddx[1];
                    data.Lww.noalias() = period_ * data.ddli_ddw[1];
                    data.Lxw.noalias() = period_ * data.ddli_dxdw[1];
                }
                else if (PS::NStages == 3)
                {
                    data.Fx.noalias() =
                        period_ / typename PS::NumScalar(4.) * (data.dki_dx[0] + typename PS::NumScalar(3.) * data.dki_dx[2]);
                    data.Fw.noalias() =
                        period_ / typename PS::NumScalar(4.) * (data.dki_dw[0] + typename PS::NumScalar(3.) * data.dki_dw[2]);
                    data.Lx.noalias() =
                        period_ / typename PS::NumScalar(4.) * (data.dli_dx[0] + typename PS::NumScalar(3.) * data.dli_dx[2]);
                    data.Lw.noalias() =
                        period_ / typename PS::NumScalar(4.) * (data.dli_dw[0] + typename PS::NumScalar(3.) * data.dli_dw[2]);
                    data.Lxx.noalias() = period_ / typename PS::NumScalar(4.) *
                                         (data.ddli_ddx[0] + typename PS::NumScalar(3.) * data.ddli_ddx[2]);
                    data.Lww.noalias() = period_ / typename PS::NumScalar(4.) *
                                         (data.ddli_ddw[0] + typename PS::NumScalar(3.) * data.ddli_ddw[2]);
                    data.Lxw.noalias() = period_ / typename PS::NumScalar(4.) *
                                         (data.ddli_dxdw[0] + typename PS::NumScalar(3.) * data.ddli_dxdw[2]);
                }
                else
                {
                    data.Fx.noalias() = period_ / typename PS::NumScalar(6.) *
                                        (data.dki_dx[0] + typename PS::NumScalar(2.) * data.dki_dx[1] +
                                         typename PS::NumScalar(2.) * data.dki_dx[2] + data.dki_dx[3]);
                    data.Fw.noalias() = period_ / typename PS::NumScalar(6.) *
                                        (data.dki_dw[0] + typename PS::NumScalar(2.) * data.dki_dw[1] +
                                         typename PS::NumScalar(2.) * data.dki_dw[2] + data.dki_dw[3]);
                    data.Lx.noalias() = period_ / typename PS::NumScalar(6.) *
                                        (data.dli_dx[0] + typename PS::NumScalar(2.) * data.dli_dx[1] +
                                         typename PS::NumScalar(2.) * data.dli_dx[2] + data.dli_dx[3]);
                    data.Lw.noalias() = period_ / typename PS::NumScalar(6.) *
                                        (data.dli_dw[0] + typename PS::NumScalar(2.) * data.dli_dw[1] +
                                         typename PS::NumScalar(2.) * data.dli_dw[2] + data.dli_dw[3]);
                    data.Lxx.noalias() = period_ / typename PS::NumScalar(6.) *
                                         (data.ddli_ddx[0] + typename PS::NumScalar(2.) * data.ddli_ddx[1] +
                                          typename PS::NumScalar(2.) * data.ddli_ddx[2] + data.ddli_ddx[3]);
                    data.Lww.noalias() = period_ / typename PS::NumScalar(6.) *
                                         (data.ddli_ddw[0] + typename PS::NumScalar(2.) * data.ddli_ddw[1] +
                                          typename PS::NumScalar(2.) * data.ddli_ddw[2] + data.ddli_ddw[3]);
                    data.Lxw.noalias() = period_ / typename PS::NumScalar(6.) *
                                         (data.ddli_dxdw[0] + typename PS::NumScalar(2.) * data.ddli_dxdw[1] +
                                          typename PS::NumScalar(2.) * data.ddli_dxdw[2] + data.ddli_dxdw[3]);
                }
                data.Hx = k0_data.Hx;
                data.Gx = k0_data.Gx;
                data.Hw.resize(node_.get_nh(), PS::NW);
                data.Gw.resize(node_.get_ng(), PS::NW);
                control_.multiplyByJacobian(w0_data, k0_data.Hu, data.Hw);
                control_.multiplyByJacobian(w0_data, k0_data.Gu, data.Gw);

                state_.JintegrateTransport(x, data.dx, data.Fx, second);
                state_.Jintegrate(x, data.dx, data.Fx, data.Fx, first, addto);
                state_.JintegrateTransport(x, data.dx, data.Fw, second);
            }

            template <typename StateVectorType, typename ControlParamVectorType>
            void quasiStatic(typename PS::SegmentData_t &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlParamVectorType> &ws,
                             const std::size_t maxiter, const typename PS::NumScalar tol) const
            {
                typename PS::ControlParamData_t &w0_data = data.controls[0];
                w0_data.u *= 0.;
                node_.quasiStatic(data.nodes[0], w0_data.u, x.derived(), maxiter, tol);
                control_.params(w0_data, 0., w0_data.u);
                w = w0_data.w;
            }

        protected:
            typename PS::State_t state_;
            typename PS::ControlParamModel_t control_;
            typename PS::NodeModel_t node_;

            typename PS::StageCoefficients_t stage_coefficients_;
            typename PS::Quadrature_t quadrature_;
            typename PS::Timings_t timings_;

            typename PS::NumScalar period_;

        }; // class SegmentERKModel234Tpl

        // A wrapper that collapses <PS, _RKType> into a single template <PS>.
        template <RKType _RKType>
        struct SegmentERKMeta234
        {
            template <typename PhaseSpec>
            struct Implementation
            {
                using Model = SegmentERKModel234Tpl<PhaseSpec, _RKType>;
                using Data = SegmentERKData234Tpl<PhaseSpec, _RKType>;
            };
        };

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_234_hpp__
