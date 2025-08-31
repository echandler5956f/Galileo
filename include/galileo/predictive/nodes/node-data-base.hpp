#ifndef __galileo_predictive_nodes_node_data_base_hpp__
#define __galileo_predictive_nodes_node_data_base_hpp__

#include "galileo/predictive/nodes/node-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct NodeDataBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using XAcc_t = ArenaMatrixTpl<typename PS::XAcc_t>;
        using XAccx_t = ArenaMatrixTpl<typename PS::XAccx_t>;
        using XAccu_t = ArenaMatrixTpl<typename PS::XAccu_t>;
        using L_t = typename PS::L_t;
        using Lx_t = ArenaMatrixTpl<typename PS::Lx_t>;
        using Lu_t = ArenaMatrixTpl<typename PS::Lu_t>;
        using Lxx_t = ArenaMatrixTpl<typename PS::Lxx_t>;
        using Lxu_t = ArenaMatrixTpl<typename PS::Lxu_t>;
        using Luu_t = ArenaMatrixTpl<typename PS::Luu_t>;
        using H_t = ArenaMatrixTpl<typename PS::H_t>;
        using Hx_t = ArenaMatrixTpl<typename PS::Hx_t>;
        using Hu_t = ArenaMatrixTpl<typename PS::Hu_t>;
        using G_t = ArenaMatrixTpl<typename PS::G_t>;
        using Gx_t = ArenaMatrixTpl<typename PS::Gx_t>;
        using Gu_t = ArenaMatrixTpl<typename PS::Gu_t>;

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
