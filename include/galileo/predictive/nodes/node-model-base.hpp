#ifndef __galileo_predictive_nodes_node_model_base_hpp__
#define __galileo_predictive_nodes_node_model_base_hpp__

#include "galileo/predictive/nodes/node-base.hpp"

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

            template <typename StateVectorType, typename ControlVectorType>
            void calc(typename PS::NodeData_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(typename PS::NodeData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void quasiStatic(typename PS::NodeData_t &data, const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlVectorType> &u,
                             const std::size_t maxiter, const typename PS::NumScalar tol) const
            {
                derived().quasiStatic(data, x.derived(), u.derived(), maxiter, tol);
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
