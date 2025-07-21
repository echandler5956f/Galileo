#ifndef __galileo_core_constraints_equality_constraint_unary_visitor_hpp__
#define __galileo_core_constraints_equality_constraint_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/common/meta/fusion.hpp"
#include "galileo/core/constraints/equality/constraint-base.hpp"

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
                typename PhaseSpec,
                template <typename PS> class ConstraintCollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
                ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>, ArgsTmp>
                    visitor(constraint_data, args);
                return boost::apply_visitor(visitor, constraint_model);
            }

            template <typename PhaseSpec, template <typename PS> class ConstraintCollectionTpl>
            static ReturnType run(
                const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
                ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
            {
                InternalVisitorModelAndData<ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>, NoArg>
                    visitor(constraint_data);
                return boost::apply_visitor(visitor, constraint_model);
            }

            template <typename ConstraintModelType, typename ArgsTmp>
            static ReturnType run(
                const ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS> &constraint_model,
                typename ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS>::Data_t &constraint_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ConstraintModelType, ArgsTmp> visitor(constraint_data, args);
                return visitor(constraint_model.derived());
            }

            template <typename ConstraintModelType>
            static ReturnType run(
                const ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS> &constraint_model,
                typename ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS>::Data_t &constraint_data)
            {
                InternalVisitorModelAndData<ConstraintModelType, NoArg> visitor(constraint_data);
                return visitor(constraint_model.derived());
            }

            template <
                typename PhaseSpec,
                template <typename PS> class ConstraintCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, constraint_model);
            }

            template <
                typename PhaseSpec,
                template <typename PS> class ConstraintCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, constraint_data);
            }

            template <typename PhaseSpec, template <typename PS> class ConstraintCollectionTpl>
            static ReturnType run(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, constraint_model);
            }

            template <typename PhaseSpec, template <typename PS> class ConstraintCollectionTpl>
            static ReturnType run(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, constraint_data);
            }

            template <typename ConstraintModelType, typename ArgsTmp>
            static ReturnType run(const ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS> &constraint_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(constraint_model.derived());
            }

            template <typename ConstraintDataType, typename ArgsTmp>
            static ReturnType run(const ConstraintDataBase<ConstraintDataType, typename traits<ConstraintDataType>::PS> &constraint_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(constraint_data.derived());
            }

            template <typename ConstraintModelType>
            static ReturnType run(const ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS> &constraint_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(constraint_model.derived());
            }

            template <typename ConstraintDataType>
            static ReturnType run(const ConstraintDataBase<ConstraintDataType, typename traits<ConstraintDataType>::PS> &constraint_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(constraint_data.derived());
            }

        private:
            template <typename ConstraintModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using ConstraintData = typename traits<ConstraintModel>::Data_t;

                InternalVisitorModelAndData(ConstraintData &constraint_data, ArgType args)
                    : constraint_data(constraint_data), args(args)
                {
                }

                template <typename ConstraintModelType>
                ReturnType operator()(const ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS> &constraint_model) const
                {
                    return bf::invoke(
                        &ConstraintVisitorDerived::template algo<ConstraintModelType>,
                        gf::append(
                            boost::ref(constraint_model.derived()),
                            boost::ref(
                                boost::get<typename ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS>::Data_t>(constraint_data)),
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
                using ConstraintData = typename traits<ConstraintModel>::Data_t;

                InternalVisitorModelAndData(ConstraintData &constraint_data)
                    : constraint_data(constraint_data)
                {
                }

                template <typename ConstraintModelType>
                ReturnType operator()(const ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS> &constraint_model) const
                {
                    return bf::invoke(
                        &ConstraintVisitorDerived::template algo<ConstraintModelType>,
                        bf::make_vector(
                            boost::ref(constraint_model.derived()),
                            boost::ref(
                                boost::get<typename ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS>::Data_t>(constraint_data))));
                }

                ReturnType operator()(const ConstraintModelVoid)
                {
                    return;
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

                template <typename ConstraintModelType>
                ReturnType operator()(const ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS> &constraint_model) const
                {
                    return bf::invoke(
                        &ConstraintVisitorDerived::template algo<ConstraintModelType>,
                        gf::append(boost::ref(constraint_model.derived()), args));
                }

                template <typename ConstraintDataType>
                ReturnType operator()(const ConstraintDataBase<ConstraintDataType, typename traits<ConstraintDataType>::PS> &constraint_data) const
                {
                    return bf::invoke(
                        &ConstraintVisitorDerived::template algo<ConstraintDataType>,
                        gf::append(boost::ref(constraint_data.derived()), args));
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

                template <typename ConstraintModelType>
                ReturnType operator()(const ConstraintModelBase<ConstraintModelType, typename traits<ConstraintModelType>::PS> &constraint_model) const
                {
                    return ConstraintVisitorDerived::template algo<ConstraintModelType>(constraint_model.derived());
                }

                template <typename ConstraintDataType>
                ReturnType operator()(const ConstraintDataBase<ConstraintDataType, typename traits<ConstraintDataType>::PS> &constraint_data) const
                {
                    return ConstraintVisitorDerived::template algo<ConstraintDataType>(constraint_data.derived());
                }
            };

        }; // struct ConstraintUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_unary_visitor_hpp__
