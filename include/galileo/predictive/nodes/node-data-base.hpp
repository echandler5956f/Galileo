#ifndef __galileo_predictive_nodes_node_data_base_hpp__
#define __galileo_predictive_nodes_node_data_base_hpp__

#include "galileo/predictive/nodes/node-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct NodeDataBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        FORWARD_ACCESSOR(CostDataManager_t, costs);
        FORWARD_ACCESSOR(ConstraintDataManager_t, constraints);
        FORWARD_ACCESSOR(XAcc_t, XAcc);
        FORWARD_ACCESSOR(XAccx_t, XAccx);
        FORWARD_ACCESSOR(XAccu_t, XAccu);
        FORWARD_ACCESSOR(L_t, L);
        FORWARD_ACCESSOR(Lx_t, Lx);
        FORWARD_ACCESSOR(Lu_t, Lu);
        FORWARD_ACCESSOR(Lxx_t, Lxx);
        FORWARD_ACCESSOR(Lxu_t, Lxu);
        FORWARD_ACCESSOR(Luu_t, Luu);
        FORWARD_ACCESSOR(H_t, H);
        FORWARD_ACCESSOR(Hx_t, Hx);
        FORWARD_ACCESSOR(Hu_t, Hu);
        FORWARD_ACCESSOR(G_t, G);
        FORWARD_ACCESSOR(Gx_t, Gx);
        FORWARD_ACCESSOR(Gu_t, Gu);

    protected:
        inline NodeDataBase() {}
        inline NodeDataBase(const NodeDataBase &clone) { *this = clone; }
        inline NodeDataBase &operator=(const NodeDataBase &clone) { return *this; }

    }; // struct NodeDataBase

} // namespace galileo

#endif // __galileo_predictive_nodes_node_data_base_hpp__
