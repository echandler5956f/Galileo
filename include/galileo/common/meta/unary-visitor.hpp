#ifndef __galileo_common_meta_unary_visitor_hpp__
#define __galileo_common_meta_unary_visitor_hpp__

#define BOOST_FUSION_INVOKE_MAX_ARITY 12

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/include/invoke.hpp>
#include <boost/variant/static_visitor.hpp>

#include <type_traits>

namespace galileo
{
    namespace fusion
    {

        namespace bf = boost::fusion;
        namespace gf = galileo::fusion; // preventing clashes with pinocchio

        using NoArg = boost::blank;

        // Append the element T at the front of boost fusion vector V.
        template <typename T, typename V>
        typename bf::result_of::push_front<V const, T>::type append(T const &t, V const &v)
        {
            return bf::push_front(v, t);
        }

        // Append the elements T1 and T2 at the front of boost fusion vector V.
        template <typename T1, typename T2, typename V>
        typename bf::result_of::push_front<typename bf::result_of::push_front<V const, T2>::type const, T1>::
            type
            append(T1 const &t1, T2 const &t2, V const &v)
        {
            return bf::push_front(bf::push_front(v, t2), t1);
        }

        // Append the elements T1, T2 and T3 at the front of boost fusion vector V.
        template <typename T1, typename T2, typename T3, typename V>
        typename bf::result_of::push_front<
            typename bf::result_of::push_front<typename bf::result_of::push_front<V const, T3>::type const, T2>::
                type const,
            T1>::type
        append(T1 const &t1, T2 const &t2, T3 const &t3, V const &v)
        {
            return bf::push_front(bf::push_front(bf::push_front(v, t3), t2), t1);
        }

        // Append the elements T1, T2, T3 and T4 at the front of boost fusion vector V.
        template <typename T1, typename T2, typename T3, typename T4, typename V>
        typename bf::result_of::push_front<
            typename bf::result_of::push_front<
                typename bf::result_of::push_front<typename bf::result_of::push_front<V const, T4>::type const, T3>::type const,
                T2>::type const,
            T1>::type
        append(T1 const &t1, T2 const &t2, T3 const &t3, T4 const &t4, V const &v)
        {
            return bf::push_front(bf::push_front(bf::push_front(bf::push_front(v, t4), t3), t2), t1);
        }

        // Append the elements T1, T2, T3, T4 and T5 at the front of boost fusion vector V.
        template <typename T1, typename T2, typename T3, typename T4, typename T5, typename V>
        typename bf::result_of::push_front<
            typename bf::result_of::push_front<
                typename bf::result_of::push_front<
                    typename bf::result_of::push_front<typename bf::result_of::push_front<V const, T5>::type const, T4>::type const,
                    T3>::type const,
                T2>::type const,
            T1>::type
        append(T1 const &t1, T2 const &t2, T3 const &t3, T4 const &t4, T5 const &t5, V const &v)
        {
            return bf::push_front(bf::push_front(bf::push_front(bf::push_front(bf::push_front(v, t5), t4), t3), t2), t1);
        }

        // Generic family traits that must be specialized for each family
        template <typename FamilyTag>
        struct UnaryVisitorFamilyTraits
        {
            template <typename PS, template <typename> class CollectionTpl>
            using ModelTpl = void; // Must be specialized

            template <typename PS, template <typename> class CollectionTpl>
            using DataTpl = void; // Must be specialized

            template <typename ModelType, typename PS>
            using ModelBase = void; // Must be specialized

            template <typename DataType, typename PS>
            using DataBase = void; // Must be specialized
        };

        // Generic unary visitor base implementation
        template <typename FamilyTag, typename VisitorDerived, typename ReturnType = void>
        struct UnaryVisitorBase
        {
        private:
            using FamilyTraits = UnaryVisitorFamilyTraits<FamilyTag>;

            template <typename PS, template <typename> class CollectionTpl>
            using ModelTpl = typename FamilyTraits::template ModelTpl<PS, CollectionTpl>;

            template <typename PS, template <typename> class CollectionTpl>
            using DataTpl = typename FamilyTraits::template DataTpl<PS, CollectionTpl>;

            template <typename ModelType, typename PS>
            using ModelBase = typename FamilyTraits::template ModelBase<ModelType, PS>;

            template <typename DataType, typename PS>
            using DataBase = typename FamilyTraits::template DataBase<DataType, PS>;

            template <typename ModelType>
            using ModelBaseOf_t = ModelBase<ModelType, typename traits<ModelType>::PS>;

            template <typename DataType>
            using DataBaseOf_t = DataBase<DataType, typename traits<DataType>::PS>;

            template <typename ModelType>
            using DataOf_t = typename ModelBaseOf_t<ModelType>::Data_t;

        public:
            // Model + Data + Args
            template <
                typename PhaseSpec,
                template <typename PS> class CollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const ModelTpl<PhaseSpec, CollectionTpl> &model,
                DataTpl<PhaseSpec, CollectionTpl> &data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ModelTpl<PhaseSpec, CollectionTpl>, ArgsTmp>
                    visitor(data, args);
                return boost::apply_visitor(visitor, model);
            }

            // Model + Data (no args)
            template <typename PhaseSpec, template <typename PS> class CollectionTpl>
            static ReturnType run(
                const ModelTpl<PhaseSpec, CollectionTpl> &model,
                DataTpl<PhaseSpec, CollectionTpl> &data)
            {
                InternalVisitorModelAndData<ModelTpl<PhaseSpec, CollectionTpl>, NoArg>
                    visitor(data);
                return boost::apply_visitor(visitor, model);
            }

            // Base Model + Data + Args
            template <typename ModelType, typename ArgsTmp>
            static ReturnType run(
                const ModelBaseOf_t<ModelType> &model,
                DataOf_t<ModelType> &data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ModelType, ArgsTmp> visitor(data, args);
                return visitor(model.derived());
            }

            // Base Model + Data (no args)
            template <typename ModelType>
            static ReturnType run(
                const ModelBaseOf_t<ModelType> &model,
                DataOf_t<ModelType> &data)
            {
                InternalVisitorModelAndData<ModelType, NoArg> visitor(data);
                return visitor(model.derived());
            }

            // Model + Args
            template <
                typename PhaseSpec,
                template <typename PS> class CollectionTpl,
                typename ArgsTmp>
            static ReturnType run(const ModelTpl<PhaseSpec, CollectionTpl> &model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, model);
            }

            // Data + Args
            template <
                typename PhaseSpec,
                template <typename PS> class CollectionTpl,
                typename ArgsTmp>
            static ReturnType run(const DataTpl<PhaseSpec, CollectionTpl> &data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, data);
            }

            // Model only
            template <typename PhaseSpec, template <typename PS> class CollectionTpl>
            static ReturnType run(const ModelTpl<PhaseSpec, CollectionTpl> &model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, model);
            }

            // Data only
            template <typename PhaseSpec, template <typename PS> class CollectionTpl>
            static ReturnType run(const DataTpl<PhaseSpec, CollectionTpl> &data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, data);
            }

            // Base Model + Args
            template <typename ModelType, typename ArgsTmp>
            static ReturnType run(const ModelBaseOf_t<ModelType> &model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(model.derived());
            }

            // Base Data + Args
            template <typename DataType, typename ArgsTmp>
            static ReturnType run(const DataBaseOf_t<DataType> &data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(data.derived());
            }

            // Base Model only
            template <typename ModelType>
            static ReturnType run(const ModelBaseOf_t<ModelType> &model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(model.derived());
            }

            // Base Data only
            template <typename DataType>
            static ReturnType run(const DataBaseOf_t<DataType> &data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(data.derived());
            }

        private:
            // Internal visitor for Model + Data + Args
            template <typename Model, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using Data = typename traits<Model>::Data_t;

                InternalVisitorModelAndData(Data &data_, ArgType args_)
                    : data(data_), args(args_) {}

                template <typename ModelType>
                ReturnType operator()(const ModelBaseOf_t<ModelType> &model) const
                {
                    return bf::invoke(
                        &VisitorDerived::template algo<ModelType>,
                        gf::append(
                            boost::ref(model.derived()),
                            boost::ref(boost::get<DataOf_t<ModelType>>(data)),
                            args));
                }

                Data &data;
                ArgType args;
            };

            // Specialization for NoArg
            template <typename Model>
            struct InternalVisitorModelAndData<Model, NoArg> : public boost::static_visitor<ReturnType>
            {
                using Data = typename traits<Model>::Data_t;

                InternalVisitorModelAndData(Data &data_) : data(data_) {}

                template <typename ModelType>
                ReturnType operator()(const ModelBaseOf_t<ModelType> &model) const
                {
                    return bf::invoke(
                        &VisitorDerived::template algo<ModelType>,
                        bf::make_vector(
                            boost::ref(model.derived()),
                            boost::ref(boost::get<DataOf_t<ModelType>>(data))));
                }

                Data &data;
            };

            // Internal visitor for Model/Data only + Args
            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args_) : args(args_) {}

                template <typename ModelType>
                ReturnType operator()(const ModelBaseOf_t<ModelType> &model) const
                {
                    return bf::invoke(
                        &VisitorDerived::template algo<ModelType>,
                        gf::append(boost::ref(model.derived()), args));
                }

                template <typename DataType>
                ReturnType operator()(const DataBaseOf_t<DataType> &data) const
                {
                    return bf::invoke(
                        &VisitorDerived::template algo<DataType>,
                        gf::append(boost::ref(data.derived()), args));
                }

                ArgType args;
            };

            // Specialization for NoArg
            template <typename Dummy>
            struct InternalVisitorModel<NoArg, Dummy> : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel() {}

                template <typename ModelType>
                ReturnType operator()(const ModelBaseOf_t<ModelType> &model) const
                {
                    return VisitorDerived::template algo<ModelType>(model.derived());
                }

                template <typename DataType>
                ReturnType operator()(const DataBaseOf_t<DataType> &data) const
                {
                    return VisitorDerived::template algo<DataType>(data.derived());
                }
            };

        }; // struct UnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_common_meta_unary_visitor_hpp__
