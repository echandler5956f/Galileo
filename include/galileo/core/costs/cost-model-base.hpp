#ifndef __galileo_core_costs_cost_model_base_hpp__
#define __galileo_core_costs_cost_model_base_hpp__

#include "galileo/core/costs/cost-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class CostModelBase : public internal::CRTP<Derived>
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
            this->derived().calcDiff(data, x.derived());
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return this->derived().createData(collector);
        }

        const PS &get_ps() const
        {
            return ps_;
        }

    protected:
        inline CostModelBase(const PS &ps)
            : ps_(ps)
        {
        }

        inline CostModelBase(const CostModelBase &clone)
        {
            *this = clone;
        }

        inline CostModelBase &operator=(const CostModelBase &clone)
        {
            return *this;
        }

        const PS &ps_;

    }; // class CostModelBase

} // namespace galileo

#endif // __galileo_core_costs_cost_model_base_hpp__
