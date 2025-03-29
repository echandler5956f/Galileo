#ifndef __galileo_predictive_phases_phase_unary_visitor_hpp__
#define __galileo_predictive_phases_phase_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/utils/fusion.hpp"
#include "galileo/predictive/phases/phase-base.hpp"

namespace galileo
{
    namespace fusion
    {

        // Base structure for Unary visitation of a PhaseModel.
        // This structure provides runners to call the right visitor according to the number of
        // arguments.
        template <typename PhaseVisitorDerived, typename ReturnType = void>
        struct PhaseUnaryVisitorBase
        {
            template <
                typename VarScalar,
                typename NumScalar,
                int Options,
                template <typename, typename, int> class PhaseCollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const predictive::PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
                predictive::PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<predictive::PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>, ArgsTmp>
                    visitor(phase_data, args);
                return boost::apply_visitor(visitor, phase_model);
            }

            template <typename VarScalar, typename NumScalar, int Options, template <typename, typename, int> class PhaseCollectionTpl>
            static ReturnType run(
                const predictive::PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
                predictive::PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data)
            {
                InternalVisitorModelAndData<predictive::PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>, NoArg>
                    visitor(phase_data);
                return boost::apply_visitor(visitor, phase_model);
            }

            template <typename PhaseModelDerived, typename ArgsTmp>
            static ReturnType run(
                const predictive::PhaseModelBase<PhaseModelDerived> &phase_model,
                typename predictive::PhaseModelBase<PhaseModelDerived>::PhaseDataDerived &phase_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<PhaseModelDerived, ArgsTmp> visitor(phase_data, args);
                return visitor(phase_model.derived());
            }

            template <typename PhaseModelDerived>
            static ReturnType run(
                const predictive::PhaseModelBase<PhaseModelDerived> &phase_model,
                typename predictive::PhaseModelBase<PhaseModelDerived>::PhaseDataDerived &phase_data)
            {
                InternalVisitorModelAndData<PhaseModelDerived, NoArg> visitor(phase_data);
                return visitor(phase_model.derived());
            }

            template <
                typename VarScalar,
                typename NumScalar,
                int Options,
                template <typename, typename, int> class PhaseCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const predictive::PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, phase_model);
            }

            template <
                typename Scalar,
                int Options,
                template <typename, int> class PhaseCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const predictive::PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, phase_data);
            }

            template <typename VarScalar, typename NumScalar, int Options, template <typename, typename, int> class PhaseCollectionTpl>
            static ReturnType run(const predictive::PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, phase_model);
            }

            template <typename VarScalar, typename NumScalar, int Options, template <typename, typename, int> class PhaseCollectionTpl>
            static ReturnType run(const predictive::PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, phase_data);
            }

            template <typename PhaseModelDerived, typename ArgsTmp>
            static ReturnType run(const predictive::PhaseModelBase<PhaseModelDerived> &phase_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(phase_model.derived());
            }

            template <typename PhaseDataDerived, typename ArgsTmp>
            static ReturnType run(const predictive::PhaseDataBase<PhaseDataDerived> &phase_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(phase_data.derived());
            }

            template <typename PhaseModelDerived>
            static ReturnType run(const predictive::PhaseModelBase<PhaseModelDerived> &phase_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(phase_model.derived());
            }

            template <typename PhaseDataDerived>
            static ReturnType run(const predictive::PhaseDataBase<PhaseDataDerived> &phase_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(phase_data.derived());
            }

        private:
            template <typename PhaseModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using PhaseData = typename PhaseModel::PhaseDataDerived;

                InternalVisitorModelAndData(PhaseData &phase_data, ArgType args)
                    : phase_data(phase_data), args(args)
                {
                }

                template <typename PhaseModelDerived>
                ReturnType operator()(const predictive::PhaseModelBase<PhaseModelDerived> &phase_model) const
                {
                    return bf::invoke(
                        &PhaseVisitorDerived::template algo<PhaseModelDerived>,
                        bf::append(
                            boost::ref(phase_model.derived()),
                            boost::ref(
                                boost::get<typename predictive::PhaseModelBase<PhaseModelDerived>::PhaseDataDerived>(phase_data)),
                            args));
                }

                ReturnType operator()(const PhaseModelVoid)
                {
                    return;
                }

                PhaseData &phase_data;
                ArgType args;
            };

            template <typename PhaseModel>
            struct InternalVisitorModelAndData<PhaseModel, NoArg>
                : public boost::static_visitor<ReturnType>
            {
                using PhaseData = typename PhaseModel::PhaseDataDerived;

                InternalVisitorModelAndData(PhaseData &phase_data)
                    : phase_data(phase_data)
                {
                }

                template <typename PhaseModelDerived>
                ReturnType operator()(const predictive::PhaseModelBase<PhaseModelDerived> &phase_model) const
                {
                    return bf::invoke(
                        &PhaseVisitorDerived::template algo<PhaseModelDerived>,
                        bf::make_vector(
                            boost::ref(phase_model.derived()),
                            boost::ref(
                                boost::get<typename predictive::PhaseModelBase<PhaseModelDerived>::PhaseDataDerived>(phase_data))));
                }

                PhaseData &phase_data;
            };

            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args)
                    : args(args)
                {
                }

                template <typename PhaseModelDerived>
                ReturnType operator()(const predictive::PhaseModelBase<PhaseModelDerived> &phase_model) const
                {
                    return bf::invoke(
                        &PhaseVisitorDerived::template algo<PhaseModelDerived>,
                        bf::append(boost::ref(phase_model.derived()), args));
                }

                template <typename PhaseDataDerived>
                ReturnType operator()(const predictive::PhaseDataBase<PhaseDataDerived> &phase_data) const
                {
                    return bf::invoke(
                        &PhaseVisitorDerived::template algo<PhaseDataDerived>,
                        bf::append(boost::ref(phase_data.derived()), args));
                }

                ReturnType operator()(const PhaseModelVoid)
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

                template <typename PhaseModelDerived>
                ReturnType operator()(const predictive::PhaseModelBase<PhaseModelDerived> &phase_model) const
                {
                    return PhaseVisitorDerived::template algo<PhaseModelDerived>(phase_model.derived());
                }

                template <typename PhaseDataDerived>
                ReturnType operator()(const predictive::PhaseDataBase<PhaseDataDerived> &phase_data) const
                {
                    return PhaseVisitorDerived::template algo<PhaseDataDerived>(phase_data.derived());
                }
            };

        }; // struct PhaseUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_predictive_phases_phase_unary_visitor_hpp__