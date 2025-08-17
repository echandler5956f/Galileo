#ifndef __galileo_core_costs_cost_model_base_hpp__
#define __galileo_core_costs_cost_model_base_hpp__

#include "galileo/core/costs/cost-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class CostModelBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calc(data, x, u);
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x, u);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return this->derived().createData(collector);
        }

    protected:
        inline CostModelBase() {}
        inline CostModelBase(const CostModelBase &clone) {}
        inline CostModelBase &operator=(const CostModelBase &clone) { return *this; }

    }; // class CostModelBase

} // namespace galileo

#endif // __galileo_core_costs_cost_model_base_hpp__
