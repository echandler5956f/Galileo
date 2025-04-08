#ifndef __galileo_predictive_phases_phase_model_base_hpp__
#define __galileo_predictive_phases_phase_model_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{
    namespace predictive
    {

        template <typename Derived, typename PhaseSpec>
        class PhaseModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using PhaseDerived = typename traits<Derived>::PhaseDerived;
            using PhaseDataDerived = typename traits<PhaseDerived>::PhaseDataDerived;
            using PhaseModelDerived = typename traits<PhaseDerived>::PhaseModelDerived;

            template <typename StateMatrixType, typename ControlParamMatrixType>
            void calc(PhaseDataDerived &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
            {
                derived().calc(data, xs.derived(), ws.derived());
            }

            template <typename StateMatrixType, typename ControlParamMatrixType>
            void calcDiff(PhaseDataDerived &data,
                          const Eigen::MatrixBase<StateMatrixType> &xs,
                          const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
            {
                derived().calcDiff(data, xs.derived(), ws.derived());
            }

            template <typename StateMatrixType, typename ControlParamMatrixType>
            void quasiStatic(PhaseDataDerived &data, const Eigen::MatrixBase<StateMatrixType> &xs,
                             Eigen::MatrixBase<ControlParamMatrixType> &ws,
                             const std::size_t maxiter, const typename PS::NumScalar &tol) const
            {
                derived().quasiStatic(data, xs.derived(), ws.derived(), maxiter, tol);
            }

            const typename PS::SegmentModel_t &segment() const
            {
                return derived().segment();
            }

            const typename PS::NumScalar &period() const
            {
                return derived().period();
            }

            int NQb() const
            {
                return derived().NQb_impl();
            }

            int NQb_impl() const
            {
                return PS::NQb;
            }

            int NQj() const
            {
                return derived().NQj_impl();
            }

            int NQj_impl() const
            {
                return PS::NQj;
            }

            int NVb() const
            {
                return derived().NVb_impl();
            }

            int NVb_impl() const
            {
                return PS::NVb;
            }

            int NVj() const
            {
                return derived().NVj_impl();
            }

            int NVj_impl() const
            {
                return PS::NVj;
            }

            int NRotors() const
            {
                return derived().NRotors_impl();
            }

            int NRotors_impl() const
            {
                return PS::NRotors;
            }

            int NQ() const
            {
                return derived().NQ_impl();
            }

            int NQ_impl() const
            {
                return PS::NQ;
            }

            int NV() const
            {
                return derived().NV_impl();
            }

            int NV_impl() const
            {
                return PS::NV;
            }

            int NX() const
            {
                return derived().NX_impl();
            }

            int NX_impl() const
            {
                return PS::NX;
            }

            int NDX() const
            {
                return derived().NDX_impl();
            }

            int NDX_impl() const
            {
                return PS::NDX;
            }

            int NUa() const
            {
                return derived().NUa_impl();
            }

            int NUa_impl() const
            {
                return PS::NUa;
            }

            int NU() const
            {
                return derived().NU_impl();
            }

            int NU_impl() const
            {
                return PS::NU;
            }

            int NOrder() const
            {
                return derived().NOrder_impl();
            }

            int NOrder_impl() const
            {
                return PS::NOrder;
            }

            int NW() const
            {
                return derived().NW_impl();
            }

            int NW_impl() const
            {
                return PS::NW;
            }

            int NStages() const
            {
                return derived().NStages_impl();
            }

            int NStages_impl() const
            {
                return PS::NStages;
            }

        protected:
            inline PhaseModelBase()
            {
            }

            inline PhaseModelBase(const PhaseModelBase &clone)
            {
                *this = clone;
            }

            inline PhaseModelBase &operator=(const PhaseModelBase &clone)
            {
                return *this;
            }

        }; // class PhaseModelBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_model_base_hpp__
