#ifndef __galileo_core_visitor_segment_unary_visitor_hpp__
#define __galileo_core_visitor_segment_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/core/visitor/fusion.hpp"
#include "galileo/core/segment/segment-base.hpp"

namespace galileo
{
    namespace fusion
    {
        // Base structure for Unary visitation of a SegmentModel.
        // This structure provides runners to call the right visitor according to the number of
        // arguments.
        template <typename SegmentVisitorDerived, typename ReturnType = void>
        struct SegmentUnaryVisitorBase
        {
            template <
                typename Scalar,
                int Options,
                template <typename, int> class SegmentCollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const SegmentModelTpl<Scalar, Options, SegmentCollectionTpl> &segment_model,
                SegmentDataTpl<Scalar, Options, SegmentCollectionTpl> &segment_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<SegmentModelTpl<Scalar, Options, SegmentCollectionTpl>, ArgsTmp>
                    visitor(segment_data, args);
                return boost::apply_visitor(visitor, segment_model);
            }

            template <typename Scalar, int Options, template <typename, int> class SegmentCollectionTpl>
            static ReturnType run(
                const SegmentModelTpl<Scalar, Options, SegmentCollectionTpl> &segment_model,
                SegmentDataTpl<Scalar, Options, SegmentCollectionTpl> &segment_data)
            {
                InternalVisitorModelAndData<SegmentModelTpl<Scalar, Options, SegmentCollectionTpl>, NoArg>
                    visitor(segment_data);
                return boost::apply_visitor(visitor, segment_model);
            }

            template <typename SegmentModelDerived, typename ArgsTmp>
            static ReturnType run(
                const SegmentModelBase<SegmentModelDerived> &segment_model,
                typename SegmentModelBase<SegmentModelDerived>::SegmentDataDerived &segment_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<SegmentModelDerived, ArgsTmp> visitor(segment_data, args);
                return visitor(segment_model.derived());
            }

            template <typename SegmentModelDerived>
            static ReturnType run(
                const SegmentModelBase<SegmentModelDerived> &segment_model,
                typename SegmentModelBase<SegmentModelDerived>::SegmentDataDerived &segment_data)
            {
                InternalVisitorModelAndData<SegmentModelDerived, NoArg> visitor(segment_data);
                return visitor(segment_model.derived());
            }

            template <
                typename Scalar,
                int Options,
                template <typename, int> class SegmentCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const SegmentModelTpl<Scalar, Options, SegmentCollectionTpl> &segment_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, segment_model);
            }

            template <
                typename Scalar,
                int Options,
                template <typename, int> class SegmentCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const SegmentDataTpl<Scalar, Options, SegmentCollectionTpl> &segment_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, segment_data);
            }

            template <typename Scalar, int Options, template <typename, int> class SegmentCollectionTpl>
            static ReturnType run(const SegmentModelTpl<Scalar, Options, SegmentCollectionTpl> &segment_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, segment_model);
            }

            template <typename Scalar, int Options, template <typename, int> class SegmentCollectionTpl>
            static ReturnType run(const SegmentDataTpl<Scalar, Options, SegmentCollectionTpl> &segment_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, segment_data);
            }

            template <typename SegmentModelDerived, typename ArgsTmp>
            static ReturnType run(const SegmentModelBase<SegmentModelDerived> &segment_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(segment_model.derived());
            }

            template <typename SegmentDataDerived, typename ArgsTmp>
            static ReturnType run(const SegmentDataBase<SegmentDataDerived> &segment_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(segment_data.derived());
            }

            template <typename SegmentModelDerived>
            static ReturnType run(const SegmentModelBase<SegmentModelDerived> &segment_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(segment_model.derived());
            }

            template <typename SegmentDataDerived>
            static ReturnType run(const SegmentDataBase<SegmentDataDerived> &segment_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(segment_data.derived());
            }

        private:
            template <typename SegmentModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                typedef typename SegmentModel::SegmentDataDerived SegmentData;

                InternalVisitorModelAndData(SegmentData &segment_data, ArgType args)
                    : segment_data(segment_data), args(args)
                {
                }

                template <typename SegmentModelDerived>
                ReturnType operator()(const SegmentModelBase<SegmentModelDerived> &segment_model) const
                {
                    return bf::invoke(
                        &SegmentVisitorDerived::template algo<SegmentModelDerived>,
                        bf::append(
                            boost::ref(segment_model.derived()),
                            boost::ref(
                                boost::get<typename SegmentModelBase<SegmentModelDerived>::SegmentDataDerived>(segment_data)),
                            args));
                }

                ReturnType operator()(const SegmentModelVoid)
                {
                    return;
                }

                SegmentData &segment_data;
                ArgType args;
            };

            template <typename SegmentModel>
            struct InternalVisitorModelAndData<SegmentModel, NoArg>
                : public boost::static_visitor<ReturnType>
            {
                typedef typename SegmentModel::SegmentDataDerived SegmentData;

                InternalVisitorModelAndData(SegmentData &segment_data)
                    : segment_data(segment_data)
                {
                }

                template <typename SegmentModelDerived>
                ReturnType operator()(const SegmentModelBase<SegmentModelDerived> &segment_model) const
                {
                    return bf::invoke(
                        &SegmentVisitorDerived::template algo<SegmentModelDerived>,
                        bf::make_vector(
                            boost::ref(segment_model.derived()),
                            boost::ref(
                                boost::get<typename SegmentModelBase<SegmentModelDerived>::SegmentDataDerived>(segment_data))));
                }

                SegmentData &segment_data;
            };

            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args)
                    : args(args)
                {
                }

                template <typename SegmentModelDerived>
                ReturnType operator()(const SegmentModelBase<SegmentModelDerived> &segment_model) const
                {
                    return bf::invoke(
                        &SegmentVisitorDerived::template algo<SegmentModelDerived>,
                        bf::append(boost::ref(segment_model.derived()), args));
                }

                template <typename SegmentDataDerived>
                ReturnType operator()(const SegmentDataBase<SegmentDataDerived> &segment_data) const
                {
                    return bf::invoke(
                        &SegmentVisitorDerived::template algo<SegmentDataDerived>,
                        bf::append(boost::ref(segment_data.derived()), args));
                }

                ReturnType operator()(const SegmentModelVoid)
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

                template <typename SegmentModelDerived>
                ReturnType operator()(const SegmentModelBase<SegmentModelDerived> &segment_model) const
                {
                    return SegmentVisitorDerived::template algo<SegmentModelDerived>(segment_model.derived());
                }

                template <typename SegmentDataDerived>
                ReturnType operator()(const SegmentDataBase<SegmentDataDerived> &segment_data) const
                {
                    return SegmentVisitorDerived::template algo<SegmentDataDerived>(segment_data.derived());
                }
            };
        }; // struct SegmentUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_core_visitor_segment_unary_visitor_hpp__