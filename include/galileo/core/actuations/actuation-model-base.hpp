#ifndef __galileo_core_actuations_actuation_model_base_hpp__
#define __galileo_core_actuations_actuation_model_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

namespace galileo
{

    template <typename Derived, typename SystemSpec>
    class ActuationModelBase : public internal::CRTP<Derived>
    {
    public:
        using SS = SystemSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using State_t = typename SS::State_t;

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

        Data_t createData(MemoryArena &arena) { return this->derived().createData(arena); }

        const SS &get_ss() const { return ss_; }
        const State_t &get_state() const { return state_; }

    protected:
        inline ActuationModelBase(const SS &ss, const State_t &state) : ss_(ss), state_(state) {}
        inline ActuationModelBase(const ActuationModelBase &clone) : ss_(clone.ss_), state_(clone.state_) {}
        inline ActuationModelBase &operator=(const ActuationModelBase &clone)
        {
            ss_ = clone.ss_;
            state_ = clone.state_;
            return *this;
        }

        SS ss_;
        State_t state_;

    }; // class ActuationModelBase

} // namespace galileo

#endif // __galileo_core_actuations_actuation_model_base_hpp__
