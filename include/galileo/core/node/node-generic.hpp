#ifndef __galileo_core_node_generic_hpp__
#define __galileo_core_node_generic_hpp__

#include "galileo/core/node/node-collection.hpp"
#include "galileo/core/node/node-basic-visitors.hxx"
#include "galileo/utils/aligned-vector.hpp"

#include <boost/mpl/contains.hpp>

namespace galileo
{

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class NodeCollectionTpl = NodeCollectionDefaultTpl>
    struct NodeTpl;
    typedef NodeTpl<context::Scalar> Node;

    template <typename _Scalar, int _Options, template <typename S, int O> class NodeCollectionTpl>
    struct traits<NodeTpl<_Scalar, _Options, NodeCollectionTpl>>
    {
        enum
        {
            Options = _Options,
            NX = Eigen::Dynamic, // Dynamic because unknown at compile time
            NU = Eigen::Dynamic,
            NDX = Eigen::Dynamic,
            NH = Eigen::Dynamic,
            NG = Eigen::Dynamic
        };

        typedef _Scalar Scalar;
        typedef NodeCollectionTpl<Scalar, Options> NodeCollection;
        typedef NodeDataTpl<Scalar, Options, NodeCollectionTpl> NodeDataDerived;
        typedef NodeModelTpl<Scalar, Options, NodeCollectionTpl> NodeModelDerived;

        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> X_t;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> U_t;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> dX_t;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> H_t;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> G_t;

        typedef const dX_t &dXTypeConstRef;
        typedef dX_t &dXTypeRef;
        typedef const H_t &HTypeConstRef;
        typedef H_t &HTypeRef;
        typedef const G_t &GTypeConstRef;
        typedef G_t &GTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &FxTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &FxTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &FuTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &FuTypeRef;
        typedef const Scalar &LTypeConstRef;
        typedef Scalar &LTypeRef;
        typedef const Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &LxTypeConstRef;
        typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &LxTypeRef;
        typedef const Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &LuTypeConstRef;
        typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic> &LuTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &LxxTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &LxxTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &LxuTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &LxuTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &LuuTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &LuuTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &HxTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &HxTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &HuTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &HuTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &GxTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &GxTypeRef;
        typedef const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &GuTypeConstRef;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &GuTypeRef;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class NodeCollectionTpl>
    struct traits<NodeDataTpl<_Scalar, _Options, NodeCollectionTpl>>
    {
        typedef NodeTpl<_Scalar, _Options, NodeCollectionTpl> NodeDerived;
        typedef typename traits<NodeDerived>::Scalar Scalar;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class NodeCollectionTpl>
    struct traits<NodeModelTpl<_Scalar, _Options, NodeCollectionTpl>>
    {
        typedef NodeTpl<_Scalar, _Options, NodeCollectionTpl> NodeDerived;
        typedef typename traits<NodeDerived>::Scalar Scalar;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class NodeCollectionTpl>
    struct NodeDataTpl
        : public NodeDataBase<NodeDataTpl<_Scalar, _Options, NodeCollectionTpl>>,
          NodeCollectionTpl<_Scalar, _Options>::NodeDataVariant
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef NodeTpl<_Scalar, _Options, NodeCollectionTpl> NodeDerived;
        typedef NodeDataBase<NodeDataTpl> Base;

        GALILEO_NODE_DATA_TYPEDEF_TEMPLATE(NodeDerived);

        typedef NodeCollectionTpl<_Scalar, _Options> NodeCollection;
        typedef typename NodeCollection::NodeDataVariant NodeDataVariant;

        using Base::operator==;
        using Base::operator!=;

        NodeDataVariant &toVariant()
        {
            return *static_cast<NodeDataVariant *>(this);
        }
        const NodeDataVariant &toVariant() const
        {
            return *static_cast<const NodeDataVariant *>(this);
        }

        dX_t XDot() const
        {
            return galileo::node_xdot(*this);
        }

        Fx_t Fx() const
        {
            return galileo::node_fx(*this);
        }

        Fu_t Fu() const
        {
            return galileo::node_fu(*this);
        }

        L_t L() const
        {
            return galileo::node_l(*this);
        }

        Lx_t Lx() const
        {
            return galileo::node_lx(*this);
        }

        Lu_t Lu() const
        {
            return galileo::node_lu(*this);
        }

        Lxx_t Lxx() const
        {
            return galileo::node_lxx(*this);
        }

        Lxu_t Lxu() const
        {
            return galileo::node_lxu(*this);
        }

        Luu_t Luu() const
        {
            return galileo::node_luu(*this);
        }

        H_t H() const
        {
            return galileo::node_h(*this);
        }

        Hx_t Hx() const
        {
            return galileo::node_hx(*this);
        }

        Hu_t Hu() const
        {
            return galileo::node_hu(*this);
        }

        G_t G() const
        {
            return galileo::node_g(*this);
        }

        Gx_t Gx() const
        {
            return galileo::node_gx(*this);
        }

        Gu_t Gu() const
        {
            return galileo::node_gu(*this);
        }

        NodeDataTpl()
            : NodeDataVariant()
        {
        }

        NodeDataTpl(const NodeDataVariant &jdata_variant)
            : NodeDataVariant(jdata_variant)
        {
        }

        template <typename NodeDataDerived>
        NodeDataTpl(const NodeDataBase<NodeDataDerived> &node_data)
            : NodeCollection::NodeDataVariant((NodeDataVariant)node_data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename NodeDataVariant::types, NodeDataDerived>));
        }

        // Define all the standard accessors
        dX_t XDot_accessor() const
        {
            return Xdot();
        }

        Fx_t Fx_accessor() const
        {
            return Fx();
        }

        Fu_t Fu_accessor() const
        {
            return Fu();
        }

        L_t L_accessor() const
        {
            return L();
        }

        Lx_t Lx_accessor() const
        {
            return Lx();
        }

        Lu_t Lu_accessor() const
        {
            return Lu();
        }

        Lxx_t Lxx_accessor() const
        {
            return Lxx();
        }

        Lxu_t Lxu_accessor() const
        {
            return Lxu();
        }

        Luu_t Luu_accessor() const
        {
            return Luu();
        }

        H_t H_accessor() const
        {
            return H();
        }

        Hx_t Hx_accessor() const
        {
            return Hx();
        }

        Hu_t Hu_accessor() const
        {
            return Hu();
        }

        G_t G_accessor() const
        {
            return G();
        }

        Gx_t Gx_accessor() const
        {
            return Gx();
        }

        Gu_t Gu_accessor() const
        {
            return Gu();
        }
    };

    template <
        typename NewScalar,
        typename Scalar,
        int Options,
        template <typename S, int O> class NodeCollectionTpl>
    struct CastType<NewScalar, NodeModelTpl<Scalar, Options, NodeCollectionTpl>>
    {
        typedef NodeModelTpl<NewScalar, Options, NodeCollectionTpl> type;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class NodeCollectionTpl>
    struct NodeModelTpl
        : NodeModelBase<NodeModelTpl<_Scalar, _Options, NodeCollectionTpl>>,
          NodeCollectionTpl<_Scalar, _Options>::NodeModelVariant
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef NodeTpl<_Scalar, _Options, NodeCollectionTpl> NodeDerived;

        GALILEO_NODE_TYPEDEF_TEMPLATE(NodeDerived);

        typedef NodeCollectionTpl<Scalar, Options> NodeCollection;
        typedef typename NodeCollection::NodeDataVariant NodeDataVariant;
        typedef typename NodeCollection::NodeModelVariant NodeModelVariant;

        NodeModelTpl()
            : NodeModelVariant()
        {
        }

        NodeModelTpl(const NodeModelVariant &jmodel_variant)
            : NodeCollection::NodeModelVariant(jmodel_variant)
        {
        }

        template <typename NodeModelDerived>
        NodeModelTpl(const NodeModelBase<NodeModelDerived> &node_model)
            : NodeModelVariant((NodeModelVariant)node_model.derived())
        {
            BOOST_MPL_ASSERT(
                (boost::mpl::contains<typename NodeModelVariant::types, NodeModelDerived>));
        }

        NodeModelVariant &toVariant()
        {
            return *static_cast<NodeModelVariant *>(this);
        }

        const NodeModelVariant &toVariant() const
        {
            return *static_cast<const NodeModelVariant *>(this);
        }

        NodeDataDerived createData() const
        {
            return ::galileo::createData<Scalar, Options, NodeCollectionTpl>(*this);
        }

        template<typename StateVectorType, typename ControlVectorType>
        void calc(NodeDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &xs,
                  const Eigen::MatrixBase<ControlVectorType> &us) const
        {
            galileo::calc<Scalar, Options, NodeCollectionTpl>(*this, data, xs.derived(), us.derived());
        }

        template<typename StateVectorType, typename ControlVectorType>
        void calcDiff(NodeDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &xs,
                      const Eigen::MatrixBase<ControlVectorType> &us) const
        {
            galileo::calcDiff<Scalar, Options, NodeCollectionTpl>(*this, data, xs.derived(), us.derived());
        }

        const X_t &get_x_lb() const
        {
            return galileo::node_get_x_lb(*this);
        }

        const X_t &get_x_ub() const
        {
            return galileo::node_get_x_ub(*this);
        }

        const U_t &get_u_lb() const
        {
            return galileo::node_get_u_lb(*this);
        }

        const U_t &get_u_ub() const
        {
            return galileo::node_get_u_ub(*this);
        }

        const H_t &get_h_lb() const
        {
            return galileo::node_get_h_lb(*this);
        }

        const H_t &get_h_ub() const
        {
            return galileo::node_get_h_ub(*this);
        }

        const G_t &get_g_lb() const
        {
            return galileo::node_get_g_lb(*this);
        }

        const G_t &get_g_ub() const
        {
            return galileo::node_get_g_ub(*this);
        }

        Eigen::Index get_nh() const
        {
            return galileo::node_get_nh(*this);
        }

        Eigen::Index get_ng() const
        {
            return galileo::node_get_ng(*this);
        }

        Eigen::Index get_id() const
        {
            return galileo::node_get_id(*this);
        }

        template <typename StateBoundVectorType>
        void set_x_lb(const Eigen::MatrixBase<StateBoundVectorType> &x_lb)
        {
            galileo::node_set_x_lb(*this, x_lb);
        }

        template <typename StateBoundVectorType>
        void set_x_ub(const Eigen::MatrixBase<StateBoundVectorType> &x_ub)
        {
            galileo::node_set_x_ub(*this, x_ub);
        }

        template <typename ControlBoundVectorType>
        void set_u_lb(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb)
        {
            galileo::node_set_u_lb(*this, u_lb);
        }

        template <typename ControlBoundVectorType>
        void set_u_ub(const Eigen::MatrixBase<ControlBoundVectorType> &u_ub)
        {
            galileo::node_set_u_ub(*this, u_ub);
        }

        template <typename EqualityBoundVectorType>
        void set_h_lb(const Eigen::MatrixBase<EqualityBoundVectorType> &h_lb)
        {
            galileo::node_set_h_lb(*this, h_lb);
        }

        template <typename EqualityBoundVectorType>
        void set_h_ub(const Eigen::MatrixBase<EqualityBoundVectorType> &h_ub)
        {
            galileo::node_set_h_ub(*this, h_ub);
        }

        template <typename InequalityBoundVectorType>
        void set_g_lb(const Eigen::MatrixBase<InequalityBoundVectorType> &g_lb)
        {
            galileo::node_set_g_lb(*this, g_lb);
        }

        template <typename InequalityBoundVectorType>
        void set_g_ub(const Eigen::MatrixBase<InequalityBoundVectorType> &g_ub)
        {
            galileo::node_set_g_ub(*this, g_ub);
        }

        void set_nh(Eigen::Index nh)
        {
            galileo::node_set_nh(*this, nh);
        }

        void set_ng(Eigen::Index ng)
        {
            galileo::node_set_ng(*this, ng);
        }

        void set_id(Eigen::Index id)
        {
            galileo::node_set_id(*this, id);
        }

        /// \returns An expression of *this with the Scalar type casted to NewScalar.
        template <typename NewScalar>
        NodeModelTpl<NewScalar, Options, NodeCollectionTpl> cast() const
        {
            return cast_node<NewScalar, Scalar, Options, NodeCollectionTpl>(*this);
        }
    };

    typedef GALILEO_ALIGNED_STD_VECTOR(NodeData) NodeDataVector;
    typedef GALILEO_ALIGNED_STD_VECTOR(NodeModel) NodeModelVector;

} // namespace galileo

#endif // __galileo_core_node_generic_hpp__
