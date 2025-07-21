#ifndef __galileo_core_costs_cost_unary_visitor_hpp__
#define __galileo_core_costs_cost_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/common/meta/fusion.hpp"
#include "galileo/core/costs/cost-base.hpp"

namespace galileo
{
    namespace fusion
    {

        // Base structure for Unary visitation of a CostModel.
        // This structure provides runners to call the right visitor according to the number of
        // arguments.
        template <typename CostVisitorDerived, typename ReturnType = void>
        struct CostUnaryVisitorBase
        {
            template <
                typename PhaseSpec,
                template <typename PS> class CostCollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
                CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<CostModelTpl<PhaseSpec, CostCollectionTpl>, ArgsTmp>
                    visitor(cost_data, args);
                return boost::apply_visitor(visitor, cost_model);
            }

            template <typename PhaseSpec, template <typename PS> class CostCollectionTpl>
            static ReturnType run(
                const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
                CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
            {
                InternalVisitorModelAndData<CostModelTpl<PhaseSpec, CostCollectionTpl>, NoArg>
                    visitor(cost_data);
                return boost::apply_visitor(visitor, cost_model);
            }

            template <typename CostModelType, typename ArgsTmp>
            static ReturnType run(
                const CostModelBase<CostModelType, typename traits<CostModelType>::PS> &cost_model,
                typename CostModelBase<CostModelType, typename traits<CostModelType>::PS>::Data_t &cost_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<CostModelType, ArgsTmp> visitor(cost_data, args);
                return visitor(cost_model.derived());
            }

            template <typename CostModelType>
            static ReturnType run(
                const CostModelBase<CostModelType, typename traits<CostModelType>::PS> &cost_model,
                typename CostModelBase<CostModelType, typename traits<CostModelType>::PS>::Data_t &cost_data)
            {
                InternalVisitorModelAndData<CostModelType, NoArg> visitor(cost_data);
                return visitor(cost_model.derived());
            }

            template <
                typename PhaseSpec,
                template <typename PS> class CostCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, cost_model);
            }

            template <
                typename PhaseSpec,
                template <typename PS> class CostCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, cost_data);
            }

            template <typename PhaseSpec, template <typename PS> class CostCollectionTpl>
            static ReturnType run(const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, cost_model);
            }

            template <typename PhaseSpec, template <typename PS> class CostCollectionTpl>
            static ReturnType run(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, cost_data);
            }

            template <typename CostModelType, typename ArgsTmp>
            static ReturnType run(const CostModelBase<CostModelType, typename traits<CostModelType>::PS> &cost_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(cost_model.derived());
            }

            template <typename CostDataType, typename ArgsTmp>
            static ReturnType run(const CostDataBase<CostDataType, typename traits<CostDataType>::PS> &cost_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(cost_data.derived());
            }

            template <typename CostModelType>
            static ReturnType run(const CostModelBase<CostModelType, typename traits<CostModelType>::PS> &cost_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(cost_model.derived());
            }

            template <typename CostDataType>
            static ReturnType run(const CostDataBase<CostDataType, typename traits<CostDataType>::PS> &cost_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(cost_data.derived());
            }

        private:
            template <typename CostModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using CostData = typename traits<CostModel>::Data_t;

                InternalVisitorModelAndData(CostData &cost_data, ArgType args)
                    : cost_data(cost_data), args(args)
                {
                }

                template <typename CostModelType>
                ReturnType operator()(const CostModelBase<CostModelType, typename traits<CostModelType>::PS> &cost_model) const
                {
                    return bf::invoke(
                        &CostVisitorDerived::template algo<CostModelType>,
                        gf::append(
                            boost::ref(cost_model.derived()),
                            boost::ref(
                                boost::get<typename CostModelBase<CostModelType, typename traits<CostModelType>::PS>::Data_t>(cost_data)),
                            args));
                }

                ReturnType operator()(const CostModelVoid)
                {
                    return;
                }

                CostData &cost_data;
                ArgType args;
            };

            template <typename CostModel>
            struct InternalVisitorModelAndData<CostModel, NoArg>
                : public boost::static_visitor<ReturnType>
            {
                using CostData = typename traits<CostModel>::Data_t;

                InternalVisitorModelAndData(CostData &cost_data)
                    : cost_data(cost_data)
                {
                }

                template <typename CostModelType>
                ReturnType operator()(const CostModelBase<CostModelType, typename traits<CostModelType>::PS> &cost_model) const
                {
                    return bf::invoke(
                        &CostVisitorDerived::template algo<CostModelType>,
                        bf::make_vector(
                            boost::ref(cost_model.derived()),
                            boost::ref(
                                boost::get<typename CostModelBase<CostModelType, typename traits<CostModelType>::PS>::Data_t>(cost_data))));
                }

                ReturnType operator()(const CostModelVoid)
                {
                    return;
                }

                CostData &cost_data;
            };

            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args)
                    : args(args)
                {
                }

                template <typename CostModelType>
                ReturnType operator()(const CostModelBase<CostModelType, typename traits<CostModelType>::PS> &cost_model) const
                {
                    return bf::invoke(
                        &CostVisitorDerived::template algo<CostModelType>,
                        gf::append(boost::ref(cost_model.derived()), args));
                }

                template <typename CostDataType>
                ReturnType operator()(const CostDataBase<CostDataType, typename traits<CostDataType>::PS> &cost_data) const
                {
                    return bf::invoke(
                        &CostVisitorDerived::template algo<CostDataType>,
                        gf::append(boost::ref(cost_data.derived()), args));
                }

                ReturnType operator()(const CostModelVoid)
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

                template <typename CostModelType>
                ReturnType operator()(const CostModelBase<CostModelType, typename traits<CostModelType>::PS> &cost_model) const
                {
                    return CostVisitorDerived::template algo<CostModelType>(cost_model.derived());
                }

                template <typename CostDataType>
                ReturnType operator()(const CostDataBase<CostDataType, typename traits<CostDataType>::PS> &cost_data) const
                {
                    return CostVisitorDerived::template algo<CostDataType>(cost_data.derived());
                }
            };

        }; // struct CostUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_core_costs_cost_unary_visitor_hpp__
