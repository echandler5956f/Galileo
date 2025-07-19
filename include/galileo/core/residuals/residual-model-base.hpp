#ifndef __galileo_core_residuals_residual_model_base_hpp__
#define __galileo_core_residuals_residual_model_base_hpp__

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

#define GALILEO_RESIDUAL_DATA_TYPEDEF(Residual)           \
    using R_t = typename traits<Residual>::R_t;           \
    using Rx_t = typename traits<Residual>::Rx_t;         \
    using Ru_t = typename traits<Residual>::Ru_t;         \
    using Arr_Rx_t = typename traits<Residual>::Arr_Rx_t; \
    using Arr_Ru_t = typename traits<Residual>::Arr_Ru_t;

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ResidualModelBase
        : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        static constexpr bool QDependent = traits<Meta_t>::QDependent;
        static constexpr bool VDependent = traits<Meta_t>::VDependent;
        static constexpr bool UDependent = traits<Meta_t>::UDependent;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

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
            this->derived().calcDiff(data, x.derived());
        }

        template <typename CostDataType, typename ActivationDataType, bool UpdateU = true>
        void calcCostDiff(CostDataType &cdata,
                          Data_t &rdata,
                          const ActivationDataType &adata) const
        {
            this->derived().calcCostDiffImpl<UpdateU>(cdata, rdata, adata);
        }

        template <typename CostDataType, typename ActivationDataType, bool UpdateU = true>
        void calcCostDiffImpl(CostDataType &cdata,
                              Data_t &rdata,
                              const ActivationDataType &adata) const
        {
            // This function computes the derivatives of the cost function based on a
            // Gauss-Newton approximation. We split the computation into two parts since it
            // is possible (and trivial) to know the optimal branch at compile time.

            calcCostDiffRxImpl(cdata, rdata, adata);

            if constexpr (UDependent && PS::DimNU_t::Value > 0 && UpdateU)
            {
                calcCostDiffRuImpl(cdata, rdata, adata);
            }
            else
            {
                if (UDependent && get_ps().get_nu() != 0 && UpdateU)
                {
                    calcCostDiffRuImpl(cdata, rdata, adata);
                }
            }
        }

        constexpr void calcCostDiffRxImpl(CostDataType &cdata,
                                          Data_t &rdata,
                                          const ActivationDataType &adata) const
        {
            if constexpr (QDependent && VDependent)
            {
                cdata.Lx.noalias() = rdata.Rx.transpose() * adata.Ar;
                rdata.Arr_Rx.noalias() = adata.Arr.diagonal().asDiagonal() * rdata.Rx;
                cdata.Lxx.noalias() = rdata.Rx.transpose() * rdata.Arr_Rx;
            }
            else if constexpr (QDependent)
            {
                Eigen::Block<Rx_t, DimNR_t::Value, PS::DimNV_t::Value, true> Rq =
                    leftCols(rdata.Rx, get_ps().get_nv_dim());
                head(cdata.Lx, get_ps().get_nv_dim()).noalias() = Rq.transpose() * adata.Ar;
                leftCols(rdata.Arr_Rx, get_ps().get_nv_dim()).noalias() =
                    adata.Arr.diagonal().asDiagonal() * Rq;
                topLeftCorner(cdata.Lxx, get_ps().get_nv_dim(), get_ps().get_nv_dim()).noalias() =
                    Rq.transpose() * leftCols(rdata.Arr_Rx, get_ps().get_nv_dim());
            }
            else if constexpr (VDependent)
            {
                Eigen::Block<Rx_t, DimNR_t::Value, PS::DimNV_t::Value, true> Rv =
                    rightCols(rdata.Rx, get_ps().get_nv_dim());
                tail(cdata.Lx, get_ps().get_nv_dim()).noalias() = Rv.transpose() * adata.Ar;
                rightCols(rdata.Arr_Rx, get_ps().get_nv_dim()).noalias() =
                    adata.Arr.diagonal().asDiagonal() * Rv;
                bottomRightCorner(cdata.Lxx, get_ps().get_nv_dim(), get_ps().get_nv_dim()).noalias() =
                    Rv.transpose() * rightCols(rdata.Arr_Rx, get_ps().get_nv_dim());
            }
        }

        constexpr void calcCostDiffRuImpl(CostDataType &cdata,
                                          Data_t &rdata,
                                          const ActivationDataType &adata) const
        {
            cdata.Lu.noalias() = rdata.Ru.transpose() * adata.Ar;
            rdata.Arr_Ru.noalias() = adata.Arr.diagonal().asDiagonal() * rdata.Ru;
            cdata.Luu.noalias() = rdata.Ru.transpose() * rdata.Arr_Ru;

            if constexpr (QDependent && VDependent)
                cdata.Lxu.noalias() = rdata.Rx.transpose() * rdata.Arr_Ru;
            else if constexpr (QDependent)
                topRows(cdata.Lxu, get_ps().get_nv_dim()).noalias() = Rq.transpose() * rdata.Arr_Ru;
            else if constexpr (VDependent)
                bottomRows(cdata.Lxu, get_ps().get_nv_dim()).noalias() =
                    Rv.transpose() * rdata.Arr_Ru;
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return this->derived().createData(collector);
        }

        const PS &get_ps() const
        {
            return ps_.get();
        }

        const DimNR_t &get_nr_dim() const
        {
            return nr_dim_;
        }

        int get_nr() const
        {
            return nr_dim_.value();
        }

        bool get_q_dependent() const
        {
            return QDependent;
        }

        bool get_v_dependent() const
        {
            return VDependent;
        }

        bool get_u_dependent() const
        {
            return UDependent;
        }

    protected:
        inline ResidualModelBase(const PS &ps, const DimNR_t &nr_dim)
            : ps_(ps), nr_dim_(nr_dim)
        {
        }

        inline ResidualModelBase(const ResidualModelBase &clone)
        {
            *this = clone;
        }

        inline ResidualModelBase &operator=(const ResidualModelBase &clone)
        {
            ps_ = clone.ps_;
            nr_dim_ = clone.nr_dim_;
            return *this;
        }

        std::reference_wrapper<const PS> ps_;
        DimNR_t nr_dim_;

    }; // class ResidualModelBase

} // namespace galileo

#endif // __galileo_core_residuals_residual_model_base_hpp__
