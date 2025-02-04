#ifndef __galileo_core_node_data_base_hpp__
#define __galileo_core_node_data_base_hpp__

#include "galileo/core/node/node-base.hpp"
#include "galileo/core/node/node-model-base.hpp"

namespace galileo
{

    template <typename Derived>
    struct NodeDataBase : CRTP<Derived>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using NodeDerived = typename traits<Derived>::NodeDerived;
        GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
        GALILEO_NODE_DATA_TYPEDEF(NodeDerived);

        FORWARD_GETTER(F);
        FORWARD_GETTER(Fx);
        FORWARD_GETTER(Fu);

        FORWARD_GETTER(L);
        FORWARD_GETTER(Lx);
        FORWARD_GETTER(Lu);
        FORWARD_GETTER(Lxx);
        FORWARD_GETTER(Lxu);
        FORWARD_GETTER(Luu);

        FORWARD_GETTER(H);
        FORWARD_GETTER(Hx);
        FORWARD_GETTER(Hu);

        FORWARD_GETTER(G);
        FORWARD_GETTER(Gx);
        FORWARD_GETTER(Gu);

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
            F_ = clone.F_;
            Fx_ = clone.Fx_;
            Fu_ = clone.Fu_;

            L_ = clone.L_;
            Lx_ = clone.Lx_;
            Lu_ = clone.Lu_;
            Lxx_ = clone.Lxx_;
            Lxu_ = clone.Lxu_;
            Luu_ = clone.Luu_;

            H_ = clone.H_;
            Hx_ = clone.Hx_;
            Hu_ = clone.Hu_;

            G_ = clone.G_;
            Gx_ = clone.Gx_;
            Gu_ = clone.Gu_;

            return *this;
        }

        F_t F_;   // Dynamics at x, u
        Fx_t Fx_; // Jacobian of the dynamics w.r.t. the state
        Fu_t Fu_; // Jacobian of the dynamics w.r.t. the control

        L_t L_;     // Cost at x, u
        Lx_t Lx_;   // Jacobian of the cost w.r.t. the state
        Lu_t Lu_;   // Jacobian of the cost w.r.t. the control
        Lxx_t Lxx_; // Hessian of the cost w.r.t. the state
        Lxu_t Lxu_; // Hessian of the cost w.r.t. the state and control
        Luu_t Luu_; // Hessian of the cost w.r.t. the control

        H_t H_;   // Equality constraint values
        Hx_t Hx_; // Jacobian of the equality constraint w.r.t. the state
        Hu_t Hu_; // Jacobian of the equality constraint w.r.t the control

        G_t G_;   // Inequality constraint values
        Gx_t Gx_; // Jacobian of the inequality constraint w.r.t. the state
        Gu_t Gu_; // Jacobian of the inequality constraint w.r.t. the control

    }; // struct NodeDataBase

} // namespace galileo

#endif // __galileo_core_node_data_base_hpp__