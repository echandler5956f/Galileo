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
    class ResidualModelBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

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

        template <typename CostDataType, typename ActivationDataType>
        void calcCostDiff(CostDataType &cdata,
                          Data_t &rdata,
                          const ActivationDataType &adata,
                          const bool update_u) const
        {
            this->derived().calcCostDiffImpl(cdata, rdata, adata, update_u);
        }

        template <typename CostDataType, typename ActivationDataType>
        void calcCostDiffImpl(CostDataType &cdata,
                              Data_t &rdata,
                              const ActivationDataType &adata,
                              const bool update_u) const
        {
            // This function computes the derivatives of the cost function based on a
            // Gauss-Newton approximation
            const bool is_ru = u_dependent() && get_nu() != 0 && update_u;
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
                Eigen::Block<Rx_t, DimNR_t::Value, PS::DimNV_t::Value, true> Rq =
                    leftCols(rdata.Rx, ps_.NV_dim);
                head(cdata.Lx, ps_.NV_dim).noalias() = Rq.transpose() * adata.Ar;
                leftCols(rdata.Arr_Rx, ps_.NV_dim).noalias() =
                    adata.Arr.diagonal().asDiagonal() * Rq;
                topLeftCorner(cdata.Lxx, ps_.NV_dim, ps_.NV_dim).noalias() =
                    Rq.transpose() * leftCols(rdata.Arr_Rx, ps_.NV_dim);
                if (is_ru)
                {
                    topRows(cdata.Lxu, ps_.NV_dim).noalias() = Rq.transpose() * rdata.Arr_Ru;
                }
            }
            else if (v_dependent())
            {
                Eigen::Block<Rx_t, DimNR_t::Value, PS::DimNV_t::Value, true> Rv =
                    rightCols(rdata.Rx, ps_.NV_dim);
                tail(cdata.Lx, ps_.NV_dim).noalias() = Rv.transpose() * adata.Ar;
                rightCols(rdata.Arr_Rx, ps_.NV_dim).noalias() =
                    adata.Arr.diagonal().asDiagonal() * Rv;
                bottomRightCorner(cdata.Lxx, ps_.NV_dim, ps_.NV_dim).noalias() =
                    Rv.transpose() * rightCols(rdata.Arr_Rx, ps_.NV_dim);
                if (is_ru)
                {
                    bottomRows(cdata.Lxu, ps_.NV_dim).noalias() =
                        Rv.transpose() * rdata.Arr_Ru;
                }
            }
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector)
        {
            return this->derived().createData(collector);
        }

        const bool q_dependent() const
        {
            return this->derived().q_dependent_impl();
        }

        const bool q_dependent_impl() const
        {
            return true;
        }

        const bool v_dependent() const
        {
            return this->derived().v_dependent_impl();
        }

        const bool v_dependent_impl() const
        {
            return true;
        }

        const bool u_dependent() const
        {
            return this->derived().u_dependent_impl();
        }

        const bool u_dependent_impl() const
        {
            return true;
        }

        const std::shared_ptr<State_t> &get_state() const
        {
            return this->derived().get_state();
        }

        const PS &get_ps() const
        {
            return ps_;
        }

        /**
         * @brief Return the dimension of the control space
         */
        const int get_nr() const
        {
            return this->derived().get_nr_impl();
        }

        const int get_nr_impl() const
        {
            if constexpr (DimNR_t::IsFixed)
            {
                return DimNR_t::Value;
            }
            else
            {
                return NR_dim_.value();
            }
        }

        const DimNR_t &NRDim() const
        {
            return NR_dim_;
        }

        /**
         * @brief Return the dimension of the control space
         */
        const int get_nu() const
        {
            return this->derived().get_nu_impl();
        }

        const int get_nu_impl() const
        {
            if constexpr (PS::DimNU_t::IsFixed)
            {
                return PS::DimNU_t::Value;
            }
            else
            {
                return ps_.NU_dim.value();
            }
        }

        const PS::DimNU_t &NUDim() const
        {
            return ps_.NU_dim;
        }

    protected:
        inline ResidualModelBase(const PS &ps, const DimNR_t &NR_dim) : ps_(ps), NR_dim_(NR_dim)
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

        const PS &ps_;
        DimNR_t NR_dim_;

    }; // class ResidualModelBase

} // namespace galileo

#endif // __galileo_core_residuals_residual_model_base_hpp__
