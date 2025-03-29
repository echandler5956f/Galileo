#ifndef __galileo_core_constraints_constraint_unary_visitor_hpp__
#define __galileo_core_constraints_constraint_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/utils/fusion.hpp"
#include "galileo/core/constraints/constraint-base.hpp"

namespace galileo
{
    namespace fusion
    {

        // Base structure for Unary visitation of a ConstraintModel.
        // This structure provides runners to call the right visitor according to the number of
        // arguments.
        template <typename ConstraintVisitorDerived, typename ReturnType = void>
        struct ConstraintUnaryVisitorBase
        {
            template <
                typename VarScalar,
                typename NumScalar,
                int Options,
                template <typename, typename, int> class ConstraintCollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const core::ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl> &constraint_model,
                core::ConstraintDataTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl> &constraint_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<core::ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl>, ArgsTmp>
                    visitor(constraint_data, args);
                return boost::apply_visitor(visitor, constraint_model);
            }

            template <typename VarScalar, typename NumScalar, int Options, template <typename, typename, int> class ConstraintCollectionTpl>
            static ReturnType run(
                const core::ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl> &constraint_model,
                core::ConstraintDataTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl> &constraint_data)
            {
                InternalVisitorModelAndData<core::ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl>, NoArg>
                    visitor(constraint_data);
                return boost::apply_visitor(visitor, constraint_model);
            }

            template <typename ConstraintModelDerived, typename ArgsTmp>
            static ReturnType run(
                const core::ConstraintModelBase<ConstraintModelDerived> &constraint_model,
                typename core::ConstraintModelBase<ConstraintModelDerived>::ConstraintDataDerived &constraint_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ConstraintModelDerived, ArgsTmp> visitor(constraint_data, args);
                return visitor(constraint_model.derived());
            }

            template <typename ConstraintModelDerived>
            static ReturnType run(
                const core::ConstraintModelBase<ConstraintModelDerived> &constraint_model,
                typename core::ConstraintModelBase<ConstraintModelDerived>::ConstraintDataDerived &constraint_data)
            {
                InternalVisitorModelAndData<ConstraintModelDerived, NoArg> visitor(constraint_data);
                return visitor(constraint_model.derived());
            }

            template <
                typename VarScalar,
                typename NumScalar,
                int Options,
                template <typename, typename, int> class ConstraintCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const core::ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl> &constraint_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, constraint_model);
            }

            template <
                typename Scalar,
                int Options,
                template <typename, int> class ConstraintCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const core::ConstraintDataTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl> &constraint_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, constraint_data);
            }

            template <typename VarScalar, typename NumScalar, int Options, template <typename, typename, int> class ConstraintCollectionTpl>
            static ReturnType run(const core::ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl> &constraint_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, constraint_model);
            }

            template <typename VarScalar, typename NumScalar, int Options, template <typename, typename, int> class ConstraintCollectionTpl>
            static ReturnType run(const core::ConstraintDataTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl> &constraint_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, constraint_data);
            }

            template <typename ConstraintModelDerived, typename ArgsTmp>
            static ReturnType run(const core::ConstraintModelBase<ConstraintModelDerived> &constraint_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(constraint_model.derived());
            }

            template <typename ConstraintDataDerived, typename ArgsTmp>
            static ReturnType run(const core::ConstraintDataBase<ConstraintDataDerived> &constraint_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(constraint_data.derived());
            }

            template <typename ConstraintModelDerived>
            static ReturnType run(const core::ConstraintModelBase<ConstraintModelDerived> &constraint_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(constraint_model.derived());
            }

            template <typename ConstraintDataDerived>
            static ReturnType run(const core::ConstraintDataBase<ConstraintDataDerived> &constraint_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(constraint_data.derived());
            }

        private:
            template <typename ConstraintModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using ConstraintData = typename ConstraintModel::ConstraintDataDerived;

                InternalVisitorModelAndData(ConstraintData &constraint_data, ArgType args)
                    : constraint_data(constraint_data), args(args)
                {
                }

                template <typename ConstraintModelDerived>
                ReturnType operator()(const core::ConstraintModelBase<ConstraintModelDerived> &constraint_model) const
                {
                    return bf::invoke(
                        &ConstraintVisitorDerived::template algo<ConstraintModelDerived>,
                        bf::append(
                            boost::ref(constraint_model.derived()),
                            boost::ref(
                                boost::get<typename core::ConstraintModelBase<ConstraintModelDerived>::ConstraintDataDerived>(constraint_data)),
                            args));
                }

                ReturnType operator()(const ConstraintModelVoid)
                {
                    return;
                }

                ConstraintData &constraint_data;
                ArgType args;
            };

            template <typename ConstraintModel>
            struct InternalVisitorModelAndData<ConstraintModel, NoArg>
                : public boost::static_visitor<ReturnType>
            {
                using ConstraintData = typename ConstraintModel::ConstraintDataDerived;

                InternalVisitorModelAndData(ConstraintData &constraint_data)
                    : constraint_data(constraint_data)
                {
                }

                template <typename ConstraintModelDerived>
                ReturnType operator()(const core::ConstraintModelBase<ConstraintModelDerived> &constraint_model) const
                {
                    return bf::invoke(
                        &ConstraintVisitorDerived::template algo<ConstraintModelDerived>,
                        bf::make_vector(
                            boost::ref(constraint_model.derived()),
                            boost::ref(
                                boost::get<typename core::ConstraintModelBase<ConstraintModelDerived>::ConstraintDataDerived>(constraint_data))));
                }

                ConstraintData &constraint_data;
            };

            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args)
                    : args(args)
                {
                }

                template <typename ConstraintModelDerived>
                ReturnType operator()(const core::ConstraintModelBase<ConstraintModelDerived> &constraint_model) const
                {
                    return bf::invoke(
                        &ConstraintVisitorDerived::template algo<ConstraintModelDerived>,
                        bf::append(boost::ref(constraint_model.derived()), args));
                }

                template <typename ConstraintDataDerived>
                ReturnType operator()(const core::ConstraintDataBase<ConstraintDataDerived> &constraint_data) const
                {
                    return bf::invoke(
                        &ConstraintVisitorDerived::template algo<ConstraintDataDerived>,
                        bf::append(boost::ref(constraint_data.derived()), args));
                }

                ReturnType operator()(const ConstraintModelVoid)
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

                template <typename ConstraintModelDerived>
                ReturnType operator()(const core::ConstraintModelBase<ConstraintModelDerived> &constraint_model) const
                {
                    return ConstraintVisitorDerived::template algo<ConstraintModelDerived>(constraint_model.derived());
                }

                template <typename ConstraintDataDerived>
                ReturnType operator()(const core::ConstraintDataBase<ConstraintDataDerived> &constraint_data) const
                {
                    return ConstraintVisitorDerived::template algo<ConstraintDataDerived>(constraint_data.derived());
                }
            };

        }; // struct ConstraintUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_core_constraints_constraint_unary_visitor_hpp__