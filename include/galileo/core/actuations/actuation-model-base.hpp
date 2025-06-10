#ifndef __galileo_core_actuations_actuation_model_base_hpp__
#define __galileo_core_actuations_actuation_model_base_hpp__

#include "galileo/core/actuations/actuation-base.hpp"

namespace galileo
{

    template <typename Derived, typename RobotSpec>
    class ActuationModelBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

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

        int nua() const
        {
            return this->derived().nua_impl();
        }

        int nua_impl() const
        {
            return RS::NUa;
        }

    protected:
        inline ActuationModelBase()
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

    }; // class ActuationModelBase

} // namespace galileo

#endif // __galileo_core_actuations_actuation_model_base_hpp__
