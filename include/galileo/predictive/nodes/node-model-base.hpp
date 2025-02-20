#ifndef __galileo_predictive_nodes_node_model_base_hpp__
#define __galileo_predictive_nodes_node_model_base_hpp__

#include "galileo/predictive/nodes/node-base.hpp"

#define GALILEO_NODE_BASIC_TYPEDEF(Node)                              \
    using Scalar = typename traits<Node>::Scalar;                     \
    using VarScalar = typename traits<Node>::VarScalar;               \
    constexpr int Options = traits<Node>::Options;                    \
    using NodeModelDerived = typename traits<Node>::NodeModelDerived; \
    using NodeDataDerived = typename traits<Node>::NodeDataDerived;

#define GALILEO_NODE_CONSTANTS(Node)       \
    constexpr int NX = traits<Node>::NX;   \
    constexpr int NU = traits<Node>::NU;   \
    constexpr int NDX = traits<Node>::NDX; \
    constexpr int NH = traits<Node>::NH;   \
    constexpr int NG = traits<Node>::NG;

#define GALILEO_NODE_MODEL_TYPEDEF(Node)                                                    \
    using State_t = typename traits<Node>::State_t;                                         \
    using ActuationModel_t = typename traits<Node>::ActuationModel_t;                       \
    using ConstraintModelCollection_t = typename traits<Node>::ConstraintModelCollection_t; \
    using CostModelCollection_t = typename traits<Node>::CostModelCollection_t;

#define GALILEO_NODE_DATA_TYPEDEF(Node)                                                   \
    using ActuationData_t = typename traits<Node>::ActuationData_t;                       \
    using ConstraintDataCollection_t = typename traits<Node>::ConstraintDataCollection_t; \
    using CostDataCollection_t = typename traits<Node>::CostDataCollection_t;             \
    using F_t = typename traits<Node>::F_t;                                               \
    using Fx_t = typename traits<Node>::Fx_t;                                             \
    using Fu_t = typename traits<Node>::Fu_t;                                             \
    using L_t = typename traits<Node>::L_t;                                               \
    using Lx_t = typename traits<Node>::Lx_t;                                             \
    using Lu_t = typename traits<Node>::Lu_t;                                             \
    using Lxx_t = typename traits<Node>::Lxx_t;                                           \
    using Lxu_t = typename traits<Node>::Lxu_t;                                           \
    using Luu_t = typename traits<Node>::Luu_t;                                           \
    using H_t = typename traits<Node>::H_t;                                               \
    using Hx_t = typename traits<Node>::Hx_t;                                             \
    using Hu_t = typename traits<Node>::Hu_t;                                             \
    using G_t = typename traits<Node>::G_t;                                               \
    using Gx_t = typename traits<Node>::Gx_t;                                             \
    using Gu_t = typename traits<Node>::Gu_t;

namespace galileo
{
    namespace predictive
    {

        template <typename Derived>
        class NodeModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using NodeDerived = typename traits<Derived>::NodeDerived;
            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_NODE_MODEL_TYPEDEF(NodeDerived);

            template <typename StateVectorType, typename ControlVectorType>
            void calc(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(NodeDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void quasiStatic(NodeDataDerived &data, const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlVectorType> &u,
                             const std::size_t maxiter, const Scalar tol) const
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

            State_t *state_;
            ActuationModel_t *actuation_;
            ConstraintModelCollection_t *constraints_;
            CostModelCollection_t *costs_;

        }; // class NodeModelBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_nodes_node_model_base_hpp__
