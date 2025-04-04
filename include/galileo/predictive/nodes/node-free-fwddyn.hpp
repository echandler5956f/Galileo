#ifndef __galileo_predictive_nodes_node_free_fwddyn_hpp__
#define __galileo_predictive_nodes_node_free_fwddyn_hpp__

#include "galileo/predictive/nodes/node-base.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename PhaseSpec>
        struct NodeFreeFwdDynTpl;

        template <typename PhaseSpec>
        struct traits<NodeFreeFwdDynTpl<PhaseSpec>>
        {
            using PS = PhaseSpec;

            using NodeDataDerived = NodeDataFreeFwdTpl<PhaseSpec>;
            using NodeModelDerived = NodeModelFreeFwdTpl<PhaseSpec>;

            // static constexpr int NU = ...

            // using RobotData_t = // TODO: add robot data
        };

        template <typename PhaseSpec>
        struct traits<NodeDataFreeFwdTpl<PhaseSpec>>
        {
            using NodeDerived = NodeFreeFwdDynTpl<PhaseSpec>;
        };

        template <typename PhaseSpec>
        struct traits<NodeModelFreeFwdTpl<PhaseSpec>>
        {
            using NodeDerived = NodeFreeFwdDynTpl<PhaseSpec>;
        };

        template <typename PhaseSpec>
        struct NodeDataFreeFwdTpl : public NodeDataBase<NodeDataFreeFwdTpl<PhaseSpec>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
            GALILEO_PHASE_SPEC_SCALARS_TYPEDEF(PS);
            GALILEO_PHASE_SPEC_CONSTANTS_TYPEDEF(PS);
            GALILEO_PHASE_SPEC_NODE_TYPES_TYPEDEF(PS);

            XAcc_t F;
            XAccx_t Fx;
            XAccu_t Fu;
            L_t L;
            Lx_t Lx;
            Lu_t Lu;
            Lxx_t Lxx;
            Lxu_t Lxu;
            Luu_t Luu;
            H_t H;
            Hx_t Hx;
            Hu_t Hu;
            G_t G;
            Gx_t Gx;
            Gu_t Gu;

        }; // class NodeDataFreeFwdTpl

        template <typename PhaseSpec>
        class NodeModelFreeFwdTpl : public NodeModelBase<NodeModelFreeFwdTpl<PhaseSpec>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
            GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

            template <typename StateVectorType, typename ControlVectorType>
            void calc(NodeData_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                // ...implementation
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(NodeData_t &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                // ...implementation
            }

            template <typename StateVectorType, typename ControlVectorType>
            void quasiStatic(NodeData_t &data, const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlVectorType> &u,
                             const std::size_t maxiter, const NumScalar tol) const
            {
                // ...implementation
            }

        }; // class NodeModelFreeFwdTpl

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_nodes_node_free_fwddyn_hpp__