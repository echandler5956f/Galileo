#ifndef __galileo_multibody_impulses_impulse_unary_visitor_hpp__
#define __galileo_multibody_impulses_impulse_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/common/meta/fusion.hpp"
#include "galileo/multibody/impulses/impulse-base.hpp"

namespace galileo
{
    namespace fusion
    {

        // Base structure for Unary visitation of a ImpulseModel.
        // This structure provides runners to call the right visitor according to the number of
        // arguments.
        template <typename ImpulseVisitorDerived, typename ReturnType = void>
        struct ImpulseUnaryVisitorBase
        {
            template <
                typename PhaseSpec,
                template <typename PS> class ImpulseCollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>, ArgsTmp>
                    visitor(impulse_data, args);
                return boost::apply_visitor(visitor, impulse_model);
            }

            template <typename PhaseSpec, template <typename PS> class ImpulseCollectionTpl>
            static ReturnType run(
                const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
            {
                InternalVisitorModelAndData<ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>, NoArg>
                    visitor(impulse_data);
                return boost::apply_visitor(visitor, impulse_model);
            }

            template <typename ImpulseModelType, typename ArgsTmp>
            static ReturnType run(
                const ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS> &impulse_model,
                typename ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS>::Data_t &impulse_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ImpulseModelType, ArgsTmp> visitor(impulse_data, args);
                return visitor(impulse_model.derived());
            }

            template <typename ImpulseModelType>
            static ReturnType run(
                const ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS> &impulse_model,
                typename ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS>::Data_t &impulse_data)
            {
                InternalVisitorModelAndData<ImpulseModelType, NoArg> visitor(impulse_data);
                return visitor(impulse_model.derived());
            }

            template <
                typename PhaseSpec,
                template <typename PS> class ImpulseCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, impulse_model);
            }

            template <
                typename PhaseSpec,
                template <typename PS> class ImpulseCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, impulse_data);
            }

            template <typename PhaseSpec, template <typename PS> class ImpulseCollectionTpl>
            static ReturnType run(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, impulse_model);
            }

            template <typename PhaseSpec, template <typename PS> class ImpulseCollectionTpl>
            static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, impulse_data);
            }

            template <typename ImpulseModelType, typename ArgsTmp>
            static ReturnType run(const ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS> &impulse_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(impulse_model.derived());
            }

            template <typename ImpulseDataType, typename ArgsTmp>
            static ReturnType run(const ImpulseDataBase<ImpulseDataType, typename traits<ImpulseDataType>::PS> &impulse_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(impulse_data.derived());
            }

            template <typename ImpulseModelType>
            static ReturnType run(const ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS> &impulse_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(impulse_model.derived());
            }

            template <typename ImpulseDataType>
            static ReturnType run(const ImpulseDataBase<ImpulseDataType, typename traits<ImpulseDataType>::PS> &impulse_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(impulse_data.derived());
            }

        private:
            template <typename ImpulseModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using ImpulseData = typename traits<ImpulseModel>::Data_t;

                InternalVisitorModelAndData(ImpulseData &impulse_data_, ArgType args_)
                    : impulse_data(impulse_data_), args(args_)
                {
                }

                template <typename ImpulseModelType>
                ReturnType operator()(const ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS> &impulse_model) const
                {
                    return bf::invoke(
                        &ImpulseVisitorDerived::template algo<ImpulseModelType>,
                        gf::append(
                            boost::ref(impulse_model.derived()),
                            boost::ref(
                                boost::get<typename ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS>::Data_t>(impulse_data)),
                            args));
                }

                ReturnType operator()(const ImpulseModelVoid)
                {
                    return;
                }

                ImpulseData &impulse_data;
                ArgType args;
            };

            template <typename ImpulseModel>
            struct InternalVisitorModelAndData<ImpulseModel, NoArg>
                : public boost::static_visitor<ReturnType>
            {
                using ImpulseData = typename traits<ImpulseModel>::Data_t;

                InternalVisitorModelAndData(ImpulseData &impulse_data_)
                    : impulse_data(impulse_data_)
                {
                }

                template <typename ImpulseModelType>
                ReturnType operator()(const ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS> &impulse_model) const
                {
                    return bf::invoke(
                        &ImpulseVisitorDerived::template algo<ImpulseModelType>,
                        bf::make_vector(
                            boost::ref(impulse_model.derived()),
                            boost::ref(
                                boost::get<typename ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS>::Data_t>(impulse_data))));
                }

                ImpulseData &impulse_data;
            };

            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args_)
                    : args(args_)
                {
                }

                template <typename ImpulseModelType>
                ReturnType operator()(const ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS> &impulse_model) const
                {
                    return bf::invoke(
                        &ImpulseVisitorDerived::template algo<ImpulseModelType>,
                        gf::append(boost::ref(impulse_model.derived()), args));
                }

                template <typename ImpulseDataType>
                ReturnType operator()(const ImpulseDataBase<ImpulseDataType, typename traits<ImpulseDataType>::PS> &impulse_data) const
                {
                    return bf::invoke(
                        &ImpulseVisitorDerived::template algo<ImpulseDataType>,
                        gf::append(boost::ref(impulse_data.derived()), args));
                }

                ReturnType operator()(const ImpulseModelVoid)
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

                template <typename ImpulseModelType>
                ReturnType operator()(const ImpulseModelBase<ImpulseModelType, typename traits<ImpulseModelType>::PS> &impulse_model) const
                {
                    return ImpulseVisitorDerived::template algo<ImpulseModelType>(impulse_model.derived());
                }

                template <typename ImpulseDataType>
                ReturnType operator()(const ImpulseDataBase<ImpulseDataType, typename traits<ImpulseDataType>::PS> &impulse_data) const
                {
                    return ImpulseVisitorDerived::template algo<ImpulseDataType>(impulse_data.derived());
                }
            };

        }; // struct ImpulseUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_multibody_impulses_impulse_unary_visitor_hpp__
