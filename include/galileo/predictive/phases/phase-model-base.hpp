#ifndef __galileo_predictive_phases_phase_model_base_hpp__
#define __galileo_predictive_phases_phase_model_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"

namespace galileo
{

    template <typename Derived, typename BasicSpec>
    class PhaseModelBase
        : public internal::CRTP<Derived>
    {
    public:
        using BS = BasicSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using NumScalar = typename BS::NumScalar;

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

        Data_t createData() const
        {
            return this->derived().createData();
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

} // namespace galileo

#endif // __galileo_predictive_phases_phase_model_base_hpp__
