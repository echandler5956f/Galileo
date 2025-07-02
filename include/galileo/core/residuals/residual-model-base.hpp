#ifndef __galileo_core_residuals_residual_model_base_hpp__
#define __galileo_core_residuals_residual_model_base_hpp__

#include "galileo/core/residuals/residual-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ResidualModelBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calc(data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x.derived(), u.derived());
        }

        template <typename CostDataType, typename ActivationDataType>
        void calcCostDiff(CostDataType &cdata,
                          Data_t &rdata,
                          const ActivationDataType &adata,
                          const bool update_u) const
        {
            // This function computes the derivatives of the cost function based on a
            // Gauss-Newton approximation
            const bool is_ru = u_dependent() && nu() != 0 && update_u;
            if (is_ru)
            {
                cdata.Lu.noalias() = rdata.Ru.transpose() * adata.Ar;
                rdata.Arr_Ru.noalias() = adata.Arr.diagonal().asDiagonal() * rdata.Ru;
                cdata.Luu.noalias() = rdata.Ru.transpose() * rdata.Arr_Ru;
            }
            if (q_dependent() && v_dependent())
            {
                cdata.Lx.noalias() = rdata.Rx.transpose() * adata.Ar;
                rdata.Arr_Rx.noalias() = adata.Arr.diagonal().asDiagonal() * rdata.Rx;
                cdata.Lxx.noalias() = rdata.Rx.transpose() * rdata.Arr_Rx;
                if (is_ru)
                {
                    cdata.Lxu.noalias() = rdata.Rx.transpose() * rdata.Arr_Ru;
                }
            }
            else if (q_dependent())
            {
                Eigen::Block<typename PS::MatrixX_t, Eigen::Dynamic, Eigen::Dynamic, true> Rq =
                    rdata.Rx.leftCols(PS::NV);
                cdata.Lx.head(PS::NV).noalias() = Rq.transpose() * adata.Ar;
                rdata.Arr_Rx.leftCols(PS::NV).noalias() =
                    adata.Arr.diagonal().asDiagonal() * Rq;
                cdata.Lxx.topLeftCorner(PS::NV, PS::NV).noalias() =
                    Rq.transpose() * rdata.Arr_Rx.leftCols(PS::NV);
                if (is_ru)
                {
                    cdata.Lxu.topRows(PS::NV).noalias() = Rq.transpose() * rdata.Arr_Ru;
                }
            }
            else if (v_dependent())
            {
                Eigen::Block<typename PS::MatrixX_t, Eigen::Dynamic, Eigen::Dynamic, true> Rv =
                    rdata.Rx.rightCols(PS::NV);
                cdata.Lx.tail(PS::NV).noalias() = Rv.transpose() * adata.Ar;
                rdata.Arr_Rx.rightCols(PS::NV).noalias() =
                    adata.Arr.diagonal().asDiagonal() * Rv;
                cdata.Lxx.bottomRightCorner(PS::NV, PS::NV).noalias() =
                    Rv.transpose() * rdata.Arr_Rx.rightCols(PS::NV);
                if (is_ru)
                {
                    cdata.Lxu.bottomRows(PS::NV).noalias() = Rv.transpose() * rdata.Arr_Ru;
                }
            }
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector)
        {
            return this->derived().createData(collector);
        }

        bool q_dependent() const
        {
            return this->derived().q_dependent_impl();
        }

        bool q_dependent_impl() const
        {
            return true;
        }

        bool v_dependent() const
        {
            return this->derived().v_dependent_impl();
        }

        bool v_dependent_impl() const
        {
            return true;
        }

        bool u_dependent() const
        {
            return this->derived().u_dependent_impl();
        }

        bool u_dependent_impl() const
        {
            return true;
        }

        int nr() const
        {
            return this->derived().nr_impl();
        }

        int nr_impl() const
        {
            return traits<Meta_t>::NR;
        }

        int nu() const
        {
            return this->derived().nu_impl();
        }

        // nu MUST be implemented in the derived class

    protected:
        inline ResidualModelBase()
        {
        }

        inline ResidualModelBase(const ResidualModelBase &clone)
        {
            *this = clone;
        }

        inline ResidualModelBase &operator=(const ResidualModelBase &clone)
        {
            return *this;
        }

    }; // class ResidualModelBase

} // namespace galileo

#endif // __galileo_core_residuals_residual_model_base_hpp__
