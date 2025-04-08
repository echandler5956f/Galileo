#ifndef __galileo_predictive_nodes_node_model_base_hpp__
#define __galileo_predictive_nodes_node_model_base_hpp__

#include "galileo/predictive/nodes/node-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{
    namespace predictive
    {

        template <typename Derived, typename PhaseSpec>
        class NodeModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using NodeDerived = typename traits<Derived>::NodeDerived;
            using NodeDataDerived = typename traits<NodeDerived>::NodeDataDerived;
            using NodeModelDerived = typename traits<NodeDerived>::NodeModelDerived;

            template <typename StateVectorType, typename ControlVectorType>
            void calc(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calc(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calc(data, x.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(NodeDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calcDiff(NodeDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calcDiff(data, x.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void quasiStatic(NodeDataDerived &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlVectorType> &u,
                             const std::size_t maxiter, const typename PS::NumScalar tol) const
            {
                derived().quasiStatic(data, x.derived(), u.derived(), maxiter, tol);
            }

            int nu() const
            {
                return derived().nu_impl();
            }

            int nu_impl() const
            {
                return traits<NodeDerived>::NU;
            }

        protected:
            inline NodeModelBase()
            {
            }

            inline NodeModelBase(const NodeModelBase &clone)
            {
                *this = clone;
            }

            inline NodeModelBase &operator=(const NodeModelBase &clone)
            {
                return *this;
            }

        }; // class NodeModelBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_nodes_node_model_base_hpp__
