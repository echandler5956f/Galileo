#ifndef __galileo_core_residuals_residual_model_base_hpp__
#define __galileo_core_residuals_residual_model_base_hpp__

#include "galileo/core/residuals/residual-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ResidualModelBase : internal::CRTP<Derived>
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

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x.derived(), u.derived());
        }

        template <typename CostDataType, typename ActivationDataType>
        void calcCostDiff(CostDataType &cdata,
                          Data_t &rdata,
                          const ActivationDataType &adata,
                          const bool update_u) const
        {
            this->derived().calcCostDiff(cdata, rdata, adata, update_u);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector)
        {
            return this->derived().createData(collector);
        }

        int nr() const
        {
            return this->derived().nr_impl();
        }

        int nr_impl() const
        {
            return traits<Meta_t>::NR;
        }

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
