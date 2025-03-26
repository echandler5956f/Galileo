#ifndef __galileo_predictive_nodes_node_data_base_hpp__
#define __galileo_predictive_nodes_node_data_base_hpp__

#include "galileo/predictive/nodes/node-base.hpp"
#include "galileo/predictive/nodes/node-model-base.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename Derived>
        struct NodeDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using NodeDerived = typename traits<Derived>::NodeDerived;
            
            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_NODE_DATA_TYPEDEF(NodeDerived);

            // F_t F;
            // Fx_t Fx;
            // Fu_t Fu;
            // L_t L;
            // Lx_t Lx;
            // Lu_t Lu;
            // Lxx_t Lxx;
            // Lxu_t Lxu;
            // Luu_t Luu;
            // H_t H;
            // Hx_t Hx;
            // Hu_t Hu;
            // G_t G;
            // Gx_t Gx;
            // Gu_t Gu;

            // DynamicsData_t dynamics;
            // ActuationData_t actuation;
            // ConstraintDataCollection_t constraints;
            // CostDataCollection_t costs;

        protected:
            inline NodeDataBase()
            {
            }

            inline NodeDataBase(const NodeDataBase &clone)
            {
                *this = clone;
            }

            inline NodeDataBase &operator=(const NodeDataBase &clone)
            {
                return *this;
            }

        }; // struct NodeDataBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_nodes_node_data_base_hpp__
