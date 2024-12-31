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
            NDX = Eigen::Dynamic,
            NU = Eigen::Dynamic,
            NH = Eigen::Dynamic,
            NG = Eigen::Dynamic
        };

        typedef _Scalar Scalar;
        typedef NodeCollectionTpl<Scalar, Options> NodeCollection;
        typedef NodeDataTpl<Scalar, Options, NodeCollectionTpl> NodeDataDerived;
        typedef NodeModelTpl<Scalar, Options, NodeCollectionTpl> NodeModelDerived;
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

        NodeDataTpl()
            : NodeDataVariant()
        {
        }

        NodeDataTpl(const NodeDataVariant &jdata_variant)
            : NodeDataVariant(jdata_variant)
        {
        }

        template <typename NodeDataDerived>
        NodeDataTpl(const NodeDataBase<NodeDataDerived> &jdata)
            : NodeCollection::NodeDataVariant((NodeDataVariant)jdata.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename NodeDataVariant::types, NodeDataDerived>));
        }

        // Define all the standard accessors
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
        NodeModelTpl(const NodeModelBase<NodeModelDerived> &jmodel)
            : NodeModelVariant((NodeModelVariant)jmodel.derived())
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

        /// \returns An expression of *this with the Scalar type casted to NewScalar.
        template <typename NewScalar>
        NodeModelTpl<NewScalar, Options, NodeCollectionTpl> cast() const
        {
            return cast_joint<NewScalar, Scalar, Options, NodeCollectionTpl>(*this);
        }
    };

    typedef GALILEO_ALIGNED_STD_VECTOR(NodeData) NodeDataVector;
    typedef GALILEO_ALIGNED_STD_VECTOR(NodeModel) NodeModelVector;

} // namespace galileo

#endif // __galileo_core_node_generic_hpp__
