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

            template <typename CostModelDerived, typename ArgsTmp>
            static ReturnType run(
                const CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS> &cost_model,
                typename CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS>::CostDataDerived &cost_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<CostModelDerived, ArgsTmp> visitor(cost_data, args);
                return visitor(cost_model.derived());
            }

            template <typename CostModelDerived>
            static ReturnType run(
                const CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS> &cost_model,
                typename CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS>::CostDataDerived &cost_data)
            {
                InternalVisitorModelAndData<CostModelDerived, NoArg> visitor(cost_data);
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

            template <typename CostModelDerived, typename ArgsTmp>
            static ReturnType run(const CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS> &cost_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(cost_model.derived());
            }

            template <typename CostDataDerived, typename ArgsTmp>
            static ReturnType run(const CostDataBase<CostDataDerived, typename traits<CostDataDerived>::PS> &cost_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(cost_data.derived());
            }

            template <typename CostModelDerived>
            static ReturnType run(const CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS> &cost_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(cost_model.derived());
            }

            template <typename CostDataDerived>
            static ReturnType run(const CostDataBase<CostDataDerived, typename traits<CostDataDerived>::PS> &cost_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(cost_data.derived());
            }

        private:
            template <typename CostModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using CostData = typename traits<CostModel>::CostDataDerived;

                InternalVisitorModelAndData(CostData &cost_data, ArgType args)
                    : cost_data(cost_data), args(args)
                {
                }

                template <typename CostModelDerived>
                ReturnType operator()(const CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS> &cost_model) const
                {
                    return bf::invoke(
                        &CostVisitorDerived::template algo<CostModelDerived>,
                        gf::append(
                            boost::ref(cost_model.derived()),
                            boost::ref(
                                boost::get<typename CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS>::CostDataDerived>(cost_data)),
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
                using CostData = typename traits<CostModel>::CostDataDerived;

                InternalVisitorModelAndData(CostData &cost_data)
                    : cost_data(cost_data)
                {
                }

                template <typename CostModelDerived>
                ReturnType operator()(const CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS> &cost_model) const
                {
                    return bf::invoke(
                        &CostVisitorDerived::template algo<CostModelDerived>,
                        bf::make_vector(
                            boost::ref(cost_model.derived()),
                            boost::ref(
                                boost::get<typename CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS>::CostDataDerived>(cost_data))));
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

                template <typename CostModelDerived>
                ReturnType operator()(const CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS> &cost_model) const
                {
                    return bf::invoke(
                        &CostVisitorDerived::template algo<CostModelDerived>,
                        gf::append(boost::ref(cost_model.derived()), args));
                }

                template <typename CostDataDerived>
                ReturnType operator()(const CostDataBase<CostDataDerived, typename traits<CostDataDerived>::PS> &cost_data) const
                {
                    return bf::invoke(
                        &CostVisitorDerived::template algo<CostDataDerived>,
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

                template <typename CostModelDerived>
                ReturnType operator()(const CostModelBase<CostModelDerived, typename traits<CostModelDerived>::PS> &cost_model) const
                {
                    return CostVisitorDerived::template algo<CostModelDerived>(cost_model.derived());
                }

                template <typename CostDataDerived>
                ReturnType operator()(const CostDataBase<CostDataDerived, typename traits<CostDataDerived>::PS> &cost_data) const
                {
                    return CostVisitorDerived::template algo<CostDataDerived>(cost_data.derived());
                }
            };

        }; // struct CostUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_core_costs_cost_unary_visitor_hpp__