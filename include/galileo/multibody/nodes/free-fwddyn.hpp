#ifndef __galileo_multibody_nodes_free_fwddyn_hpp__
#define __galileo_multibody_nodes_free_fwddyn_hpp__

#include "galileo/core/node/node-base.hpp"
#include <pinocchio/multibody/data.hpp>
#include <pinocchio/multibody/model.hpp>

namespace galileo
{

    template <
        typename NumScalar,
        typename VarScalar,
        template <typename, typename> class State,
        template <typename, typename> class Actuation,
        template <typename, typename, class> class ConstraintModelCollection,
        template <typename, typename, class> class CostModelCollection>
    struct NodeFreeFwddyn;

    template <
        typename _NumScalar,
        typename _VarScalar,
        template <typename, typename> class _State,
        template <typename, typename> class _Actuation,
        template <typename, typename, class> class _ConstraintModelCollection,
        template <typename, typename, class> class _CostModelCollection>
    struct traits<NodeFreeFwddyn<_NumScalar, _VarScalar, _State, _Actuation, _ConstraintModelCollection, _CostModelCollection>>
    {
        using NumScalar = _NumScalar;
        using VarScalar = _VarScalar;
        using State_t = _State<NumScalar, VarScalar>;
        using Actuation_t = _Actuation<Scalar, VarScalar>;
        using ConstraintModelCollection_t = _ConstraintModelCollection<NumScalar, VarScalar, State_t>;
        using CostModelCollection_t = _CostModelCollection<NumScalar, VarScalar, State_t>;

        using NodeModelDerived = NodeFreeFwddyn<NumScalar, VarScalar, State_t, Actuation_t, ConstraintModelCollection_t, CostModelCollection_t>;
        using NodeDataDerived = NodeDataFwddyn<NumScalar, VarScalar, State_t, Actuation_t, ConstraintModelCollection_t, CostModelCollection_t>;

        using X_Bounds_t = Eigen::Matrix<NumScalar, State_t::NX, 2>;
        using U_Bounds_t = Eigen::Matrix<NumScalar, Actuation_t::NU, 2>;
        using H_Bounds_t = Eigen::Matrix<NumScalar, ConstraintModelCollection_t::NH, 1>;
        using G_Bounds_t = Eigen::Matrix<NumScalar, ConstraintModelCollection_t::NG, 2>;

        using F_t = Eigen::Matrix<VarScalar, State_t::NDX, 1>;
        using Fx_t = Eigen::Matrix<VarScalar, State_t::NDX, State_t::NX>;
        using Fu_t = Eigen::Matrix<VarScalar, State_t::NDX, Actuation_t::NU>;
        using L_t = Eigen::Matrix<VarScalar, 1, 1>;
        using Lx_t = Eigen::Matrix<VarScalar, 1, State_t::NX>;
        using Lu_t = Eigen::Matrix<VarScalar, 1, Actuation_t::NX>;
        using Lxx_t = Eigen::Matrix<VarScalar, State_t::NX, State_t::NX>;
        using Lxu_t = Eigen::Matrix<VarScalar, State_t::NX, Actuation_t::NU>;
        using Luu_t = Eigen::Matrix<VarScalar, Actuation_t::NU, Actuation_t::NU>;
        using H_t = Eigen::Matrix<VarScalar, ConstraintModelCollection_t::NH, 1>;
        using Hx_t = Eigen::Matrix<VarScalar, ConstraintModelCollection_t::NH, State_t::NX>;
        using Hu_t = Eigen::Matrix<VarScalar, ConstraintModelCollection_t::NH, Actuation_t::NU>;
        using G_t = Eigen::Matrix<VarScalar, ConstraintModelCollection_t::NG, 1>;
        using Gx_t = Eigen::Matrix<VarScalar, ConstraintModelCollection_t::NG, State_t::NX>;
        using Gu_t = Eigen::Matrix<VarScalar, ConstraintModelCollection_t::NG, Actuation_t::NU>;
    };

    template <
        typename NumScalar,
        typename VarScalar,
        template <typename, typename> class State,
        template <typename, typename> class Actuation,
        template <typename, typename, class> class ConstraintModelCollection,
        template <typename, typename, class> class CostModelCollection>
    struct NodeDataFreeFwddyn; // fwd

    template <
        typename _NumScalar,
        typename _VarScalar,
        template <typename, typename> class _State,
        template <typename, typename> class _Actuation,
        template <typename, typename, class> class _ConstraintModelCollection,
        template <typename, typename, class> class _CostModelCollection>
    struct traits<NodeDataFreeFwddyn<_NumScalar, _VarScalar, _State, _Actuation, _ConstraintModelCollection, _CostModelCollection>>
    {
        using NodeDerived = NodeFreeFwddyn<_NumScalar, _VarScalar, _State, _Actuation, _ConstraintModelCollection, _CostModelCollection>;
    };

    template <
        typename NumScalar,
        typename VarScalar,
        template <typename, typename> class State,
        template <typename, typename> class Actuation,
        template <typename, typename, class> class ConstraintModelCollection,
        template <typename, typename, class> class CostModelCollection>
    struct NodeModelFreeFwddyn; // fwd

    template <
        typename _NumScalar,
        typename _VarScalar,
        template <typename, typename> class _State,
        template <typename, typename> class _Actuation,
        template <typename, typename, class> class _ConstraintModelCollection,
        template <typename, typename, class> class _CostModelCollection>
    struct traits<NodeModelFreeFwddyn<_NumScalar, _VarScalar, _State, _Actuation, _ConstraintModelCollection, _CostModelCollection>>
    {
        using NodeDerived = NodeFreeFwddyn<_NumScalar, _VarScalar, _State, _Actuation, _ConstraintModelCollection, _CostModelCollection>;
    };

    template <
        typename NumScalar,
        typename VarScalar,
        template <typename, typename> class State,
        template <typename, typename> class Actuation,
        template <typename, typename, class> class ConstraintModelCollection,
        template <typename, typename, class> class CostModelCollection>
    class NodeModelFreeFwddyn : NodeModelBase<NodeModelFreeFwddyn<NumScalar, VarScalar, State, Actuation, ConstraintModelCollection, CostModelCollection>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using NodeDerived = NodeFreeFwddyn<NumScalar, VarScalar, State, Actuation, ConstraintModelCollection, CostModelCollection>;
        GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
        GALILEO_NODE_ACTION_DEF_TYPEDEF(NodeDerived);
        GALILEO_NODE_MODEL_TYPEDEF(NodeDerived);

        pinocchio::ModelTpl<NumScalar> &pinocchio_;
        CostModelCollection_t &costs_;
        ConstraintModelCollection_t &constraints_;
        State_t &state_;
        Actuation_t &actuation_;

        template <typename StateVectorType, typename ControlVectorType>
        void calc(NodeDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            Eigen::VectorBlock<const StateVectorType, State_t::NQ> q = x.template head<State_t::NQ>();
            Eigen::VectorBlock<const StateVectorType, State_t::NV> v = x.template tail<State_t::NV>();

            actuation_->calc(d->multibody.actuation, x, u);

            d->xout = pinocchio::aba(pinocchio_, d->pinocchio, q, v,
                                     d->multibody.actuation->tau);
            pinocchio::updateGlobalPlacements(pinocchio_, d->pinocchio);

            d->multibody.joint->a = d->xout;
            d->multibody.joint->tau = u;
            costs_->calc(d->costs, x, u);
            d->cost = d->costs->cost;
            if (constraints_ != nullptr)
            {
                d->constraints->resize(this, d);
                constraints_->calc(d->constraints, x, u);
            }
        }
    };

} // namespace galileo

#endif // __galileo_multibody_nodes_free_fwddyn_hpp__
