#ifndef __galileo_predictive_nodes_node_model_base_hpp__
#define __galileo_predictive_nodes_node_model_base_hpp__

#include "galileo/predictive/nodes/node-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class NodeModelBase
        : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using NumScalar = typename PS::NumScalar;
        using State_t = typename PS::State_t;
        using RobotModel_t = typename PS::RobotModel_t;

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

        template <typename StateVectorType, typename ControlVectorType>
        void quasiStatic(Data_t &data,
                         const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlVectorType> &u,
                         const int maxiter, const NumScalar tol) const
        {
            this->derived().quasiStatic(data, x.derived(), u.derived(), maxiter, tol);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector)
        {
            return this->derived().createData(collector);
        }

        const PS &get_ps() const
        {
            return ps_.get();
        }

        const State_t &get_state() const
        {
            return get_ps().get_state();
        }

        const RobotModel_t &get_robot() const
        {
            return robot_.get();
        }

    protected:
        inline NodeModelBase(const PS &ps)
            : ps_(ps), robot_(ps.get_state().get_robot())
        {
        }

        inline NodeModelBase(const NodeModelBase &clone)
        {
            *this = clone;
        }

        inline NodeModelBase &operator=(const NodeModelBase &clone)
        {
            ps_ = clone.ps_;
            robot_ = clone.robot_;
            return *this;
        }

        std::reference_wrapper<const PS> ps_;
        std::reference_wrapper<const RobotModel_t> robot_;

    }; // class NodeModelBase

} // namespace galileo

#endif // __galileo_predictive_nodes_node_model_base_hpp__
