#ifndef __galileo_core_actuations_actuation_model_base_hpp__
#define __galileo_core_actuations_actuation_model_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"
#include "galileo/multibody/robot-spec.hpp"

#include <memory>

namespace galileo
{

    template <typename Derived, typename RobotSpec>
    class ActuationModelBase
        : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(RS);

        using Meta_t = typename RS::ActuationMeta_t;
        using Model_t = typename RS::ActuationModel_t;
        using Data_t = typename RS::ActuationData_t;

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
            // Nothing happens
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x.derived(), u.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            // Nothing happens
        }

        template <typename StateVectorType, typename TauVectorType>
        void commands(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<TauVectorType> &tau) const
        {
            this->derived().commands(data, x.derived(), tau.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void torqueTransform(Data_t &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().torqueTransform(data, x.derived(), u.derived());
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector)
        {
            return this->derived().createData(collector);
        }

        const std::shared_ptr<State_t> &get_state() const
        {
            return state_;
        }

        const int get_nua() const
        {
            return state_->get_nua();
        }

        const RS::DimNUa_t &get_nua_dim() const
        {
            return state_->get_nua_dim();
        }

    protected:
        inline ActuationModelBase(const std::shared_ptr<State_t> &state)
            : state_(state)
        {
        }

        inline ActuationModelBase(const ActuationModelBase &clone)
        {
            *this = clone;
        }

        inline ActuationModelBase &operator=(const ActuationModelBase &clone)
        {
            return *this;
        }

        const std::shared_ptr<State_t> &state_;

    }; // class ActuationModelBase

} // namespace galileo

#endif // __galileo_core_actuations_actuation_model_base_hpp__
