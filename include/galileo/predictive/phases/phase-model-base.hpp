#ifndef __galileo_predictive_phases_phase_model_base_hpp__
#define __galileo_predictive_phases_phase_model_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class PhaseModelBase
        : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using NumScalar = typename PS::NumScalar;

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs,
                  const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            this->derived().calc(data, xs, ws);
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            this->derived().calcDiff(data, xs, ws);
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void quasiStatic(Data_t &data, const Eigen::MatrixBase<StateMatrixType> &xs,
                         Eigen::MatrixBase<ControlParamMatrixType> &ws,
                         const int maxiter, const NumScalar &tol) const
        {
            this->derived().quasiStatic(data, xs, ws, maxiter, tol);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return this->derived().createData(collector);
        }

        const PS &get_ps() const
        {
            return this->derived().get_ps_impl();
        }

        const PS &get_ps_impl() const
        {
            return ps_.get();
        }

    protected:
        inline PhaseModelBase(const PS &ps)
            : ps_(ps)
        {
        }

        inline PhaseModelBase(const PhaseModelBase &clone)
            : ps_(clone.ps_)
        {
            *this = clone;
        }

        inline PhaseModelBase &operator=(const PhaseModelBase &clone)
        {
            ps_ = clone.ps_;
            return *this;
        }

        std::reference_wrapper<const PS> ps_;

    }; // class PhaseModelBase

} // namespace galileo

#endif // __galileo_predictive_phases_phase_model_base_hpp__
