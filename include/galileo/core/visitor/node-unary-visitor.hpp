#ifndef __galileo_core_visitor_node_unary_visitor_hpp__
#define __galileo_core_visitor_node_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/core/visitor/fusion.hpp"
#include "galileo/core/node/node-base.hpp"

namespace galileo
{
    namespace fusion
    {
        // Base structure for Unary visitation of a NodeModel.
        // This structure provides runners to call the right visitor according to the number of
        // arguments.
        template <typename NodeVisitorDerived, typename ReturnType = void>
        struct NodeUnaryVisitorBase
        {
            template <
                typename Scalar,
                int Options,
                template <typename, int> class NodeCollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<NodeModelTpl<Scalar, Options, NodeCollectionTpl>, ArgsTmp>
                    visitor(node_data, args);
                return boost::apply_visitor(visitor, node_model);
            }

            template <typename Scalar, int Options, template <typename, int> class NodeCollectionTpl>
            static ReturnType run(
                const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model,
                NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
            {
                InternalVisitorModelAndData<NodeModelTpl<Scalar, Options, NodeCollectionTpl>, NoArg>
                    visitor(node_data);
                return boost::apply_visitor(visitor, node_model);
            }

            template <typename NodeModelDerived, typename ArgsTmp>
            static ReturnType run(
                const NodeModelBase<NodeModelDerived> &node_model,
                typename NodeModelBase<NodeModelDerived>::NodeDataDerived &node_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<NodeModelDerived, ArgsTmp> visitor(node_data, args);
                return visitor(node_model.derived());
            }

            template <typename NodeModelDerived>
            static ReturnType run(
                const NodeModelBase<NodeModelDerived> &node_model,
                typename NodeModelBase<NodeModelDerived>::NodeDataDerived &node_data)
            {
                InternalVisitorModelAndData<NodeModelDerived, NoArg> visitor(node_data);
                return visitor(node_model.derived());
            }

            template <
                typename Scalar,
                int Options,
                template <typename, int> class NodeCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, node_model);
            }

            template <
                typename Scalar,
                int Options,
                template <typename, int> class NodeCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, node_data);
            }

            template <typename Scalar, int Options, template <typename, int> class NodeCollectionTpl>
            static ReturnType run(const NodeModelTpl<Scalar, Options, NodeCollectionTpl> &node_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, node_model);
            }

            template <typename Scalar, int Options, template <typename, int> class NodeCollectionTpl>
            static ReturnType run(const NodeDataTpl<Scalar, Options, NodeCollectionTpl> &node_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, node_data);
            }

            template <typename NodeModelDerived, typename ArgsTmp>
            static ReturnType run(const NodeModelBase<NodeModelDerived> &node_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(node_model.derived());
            }

            template <typename NodeDataDerived, typename ArgsTmp>
            static ReturnType run(const NodeDataBase<NodeDataDerived> &node_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(node_data.derived());
            }

            template <typename NodeModelDerived>
            static ReturnType run(const NodeModelBase<NodeModelDerived> &node_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(node_model.derived());
            }

            template <typename NodeDataDerived>
            static ReturnType run(const NodeDataBase<NodeDataDerived> &node_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(node_data.derived());
            }

        private:
            template <typename NodeModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                typedef typename NodeModel::NodeDataDerived NodeData;

                InternalVisitorModelAndData(NodeData &node_data, ArgType args)
                    : node_data(node_data), args(args)
                {
                }

                template <typename NodeModelDerived>
                ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
                {
                    return bf::invoke(
                        &NodeVisitorDerived::template algo<NodeModelDerived>,
                        bf::append(
                            boost::ref(node_model.derived()),
                            boost::ref(
                                boost::get<typename NodeModelBase<NodeModelDerived>::NodeDataDerived>(node_data)),
                            args));
                }

                ReturnType operator()(const NodeModelVoid)
                {
                    return;
                }

                NodeData &node_data;
                ArgType args;
            };

            template <typename NodeModel>
            struct InternalVisitorModelAndData<NodeModel, NoArg>
                : public boost::static_visitor<ReturnType>
            {
                typedef typename NodeModel::NodeDataDerived NodeData;

                InternalVisitorModelAndData(NodeData &node_data)
                    : node_data(node_data)
                {
                }

                template <typename NodeModelDerived>
                ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
                {
                    return bf::invoke(
                        &NodeVisitorDerived::template algo<NodeModelDerived>,
                        bf::make_vector(
                            boost::ref(node_model.derived()),
                            boost::ref(
                                boost::get<typename NodeModelBase<NodeModelDerived>::NodeDataDerived>(node_data))));
                }

                NodeData &node_data;
            };

            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args)
                    : args(args)
                {
                }

                template <typename NodeModelDerived>
                ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
                {
                    return bf::invoke(
                        &NodeVisitorDerived::template algo<NodeModelDerived>,
                        bf::append(boost::ref(node_model.derived()), args));
                }

                template <typename NodeDataDerived>
                ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
                {
                    return bf::invoke(
                        &NodeVisitorDerived::template algo<NodeDataDerived>,
                        bf::append(boost::ref(node_data.derived()), args));
                }

                ReturnType operator()(const NodeModelVoid)
                {
                    return;
                }

                ArgType args;
            };

            template <typename Dummy>
            struct InternalVisitorModel<NoArg, Dummy> : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel()
                {
                }

                template <typename NodeModelDerived>
                ReturnType operator()(const NodeModelBase<NodeModelDerived> &node_model) const
                {
                    return NodeVisitorDerived::template algo<NodeModelDerived>(node_model.derived());
                }

                template <typename NodeDataDerived>
                ReturnType operator()(const NodeDataBase<NodeDataDerived> &node_data) const
                {
                    return NodeVisitorDerived::template algo<NodeDataDerived>(node_data.derived());
                }
            };
        }; // struct NodeUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_core_visitor_node_unary_visitor_hpp__