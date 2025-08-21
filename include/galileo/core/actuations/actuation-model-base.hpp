#ifndef __galileo_core_actuations_actuation_model_base_hpp__
#define __galileo_core_actuations_actuation_model_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

namespace galileo
{

    template <typename Derived, typename Spec>
    class ActuationModelBase : public internal::CRTP<Derived>
    {
    public:
        using SS = Spec;

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
            // Nothing happens
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x, u);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            // Nothing happens
        }

        template <typename StateVectorType, typename TauVectorType>
        void commands(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<TauVectorType> &tau) const
        {
            this->derived().commands(data, x, tau);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void torqueTransform(Data_t &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().torqueTransform(data, x, u);
        }

        Data_t createData() { return this->derived().createData(); }

        const SS &get_spec() const { return spec_; }

    protected:
        inline ActuationModelBase(const SS &spec) : spec_(spec) {}
        inline ActuationModelBase(const ActuationModelBase &clone) : spec_(clone.spec_) {}
        inline ActuationModelBase &operator=(const ActuationModelBase &clone)
        {
            spec_ = clone.spec_;
            return *this;
        }

        SS spec_;

    }; // class ActuationModelBase

} // namespace galileo

#endif // __galileo_core_actuations_actuation_model_base_hpp__
