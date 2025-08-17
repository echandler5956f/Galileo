#ifndef __galileo_common_visitors_binary_visitor_hpp__
#define __galileo_common_visitors_binary_visitor_hpp__

#include "galileo/common/visitors/unary-visitor.hpp"

namespace galileo
{
    namespace fusion
    {

        // Must be template instantiated after the UnaryVisitorFamilyTraits for the particular family tag
        template <typename FamilyTag>
        struct BinaryVisitorFamilyTraits
        {
            template <typename LeftPS, template <typename> class LeftCollectionTpl>
            using LeftUnaryVisitorFamilyTraits = UnaryVisitorFamilyTraits<FamilyTag>;

            template <typename RightPS, template <typename> class RightCollectionTpl>
            using RightUnaryVisitorFamilyTraits = UnaryVisitorFamilyTraits<FamilyTag>;
        };

        // Generic binary visitor base implementation
        template <typename FamilyTag, typename VisitorDerived, typename ReturnType = void>
        struct BinaryVisitorBase
        {
        private:
            using FamilyTraits = BinaryVisitorFamilyTraits<FamilyTag>;

            using LeftUnaryTraits = typename FamilyTraits::LeftUnaryVisitorFamilyTraits;
            using RightUnaryTraits = typename FamilyTraits::RightUnaryVisitorFamilyTraits;

            template <typename LeftPS, template <typename> class LeftCollectionTpl>
            using LeftModelTpl = typename LeftUnaryTraits::template ModelTpl<LeftPS, LeftCollectionTpl>;
            template <typename RightPS, template <typename> class RightCollectionTpl>
            using RightModelTpl = typename RightUnaryTraits::template ModelTpl<RightPS, RightCollectionTpl>;

            template <typename LeftPS, template <typename> class LeftCollectionTpl>
            using LeftDataTpl = typename LeftUnaryTraits::template DataTpl<LeftPS, LeftCollectionTpl>;
            template <typename RightPS, template <typename> class RightCollectionTpl>
            using RightDataTpl = typename RightUnaryTraits::template DataTpl<RightPS, RightCollectionTpl>;

            template <typename LeftModelType, typename LeftPS>
            using LeftModelBase = typename LeftUnaryTraits::template ModelBase<LeftModelType, LeftPS>;
            template <typename RightModelType, typename RightPS>
            using RightModelBase = typename RightUnaryTraits::template ModelBase<RightModelType, RightPS>;

            template <typename LeftDataType, typename LeftPS>
            using LeftDataBase = typename LeftUnaryTraits::template DataBase<LeftDataType, LeftPS>;
            template <typename RightDataType, typename RightPS>
            using RightDataBase = typename RightUnaryTraits::template DataBase<RightDataType, RightPS>;

            template <typename LeftModelType>
            using LeftModelBaseOf_t = typename LeftUnaryTraits::template ModelBaseOf_t<LeftModelType>;
            template <typename RightModelType>
            using RightModelBaseOf_t = typename RightUnaryTraits::template ModelBaseOf_t<RightModelType>;

            template <typename LeftDataType>
            using LeftDataBaseOf_t = typename LeftUnaryTraits::template DataBaseOf_t<LeftDataType>;
            template <typename RightDataType>
            using RightDataBaseOf_t = typename RightUnaryTraits::template DataBaseOf_t<RightDataType>;

            template <typename LeftDataType>
            using LeftDataOf_t = typename LeftUnaryTraits::template DataOf_t<LeftDataType>;
            template <typename RightDataType>
            using RightDataOf_t = typename RightUnaryTraits::template DataOf_t<RightDataType>;

        public:
            // Model + Data + Args
            template <typename LeftPhaseSpec,
                      template <typename> class LeftCollectionTpl,
                      typename RightPhaseSpec,
                      template <typename> class RightCollectionTpl,
                      typename ArgsTmp>
            static ReturnType run(const LeftModelTpl<LeftPhaseSpec, LeftCollectionTpl> &left_model,
                                  const RightModelTpl<RightPhaseSpec, RightCollectionTpl> &right_model,
                                  LeftDataTpl<LeftPhaseSpec, LeftCollectionTpl> &left_data,
                                  RightDataTpl<RightPhaseSpec, RightCollectionTpl> &right_data,
                                  ArgsTmp args)
            {
                InternalVisitorModelAndData<LeftModelTpl<LeftPhaseSpec, LeftCollectionTpl>,
                                            RightModelTpl<RightPhaseSpec, RightCollectionTpl>,
                                            ArgsTmp>
                    visitor(left_data, right_data, args);
                return boost::apply_visitor(visitor, left_model, right_model);
            }

            // Model + Data (no args)
            template <typename LeftPhaseSpec,
                      template <typename> class LeftCollectionTpl,
                      typename RightPhaseSpec,
                      template <typename> class RightCollectionTpl>
            static ReturnType run(const LeftModelTpl<LeftPhaseSpec, LeftCollectionTpl> &left_model,
                                  const RightModelTpl<RightPhaseSpec, RightCollectionTpl> &right_model,
                                  LeftDataTpl<LeftPhaseSpec, LeftCollectionTpl> &left_data,
                                  RightDataTpl<RightPhaseSpec, RightCollectionTpl> &right_data)
            {
                InternalVisitorModelAndData<LeftModelTpl<LeftPhaseSpec, LeftCollectionTpl>,
                                            RightModelTpl<RightPhaseSpec, RightCollectionTpl>,
                                            NoArg>
                    visitor(left_data, right_data);
                return boost::apply_visitor(visitor, left_model, right_model);
            }

            // Base Model + Data + Args
            template <typename LeftModelType, typename RightModelType, typename ArgsTmp>
            static ReturnType run(const LeftModelBaseOf_t<LeftModelType> &left_model,
                                  const RightModelBaseOf_t<RightModelType> &right_model,
                                  LeftDataOf_t<LeftModelType> &left_data,
                                  RightDataOf_t<RightModelType> &right_data,
                                  ArgsTmp args)
            {
                InternalVisitorModelAndData<LeftModelType, RightModelType, ArgsTmp> visitor(
                    left_data, right_data, args);
                return visitor(left_model.derived(), right_model.derived());
            }

            // Base Model + Data (no args)
            template <typename LeftModelType, typename RightModelType>
            static ReturnType run(const LeftModelBaseOf_t<LeftModelType> &left_model,
                                  const RightModelBaseOf_t<RightModelType> &right_model,
                                  LeftDataOf_t<LeftModelType> &left_data,
                                  RightDataOf_t<RightModelType> &right_data)
            {
                InternalVisitorModelAndData<LeftModelType, RightModelType, NoArg> visitor(left_data, right_data);
                return visitor(left_model.derived(), right_model.derived());
            }

            // Model + Args
            template <typename LeftPhaseSpec,
                      template <typename> class LeftCollectionTpl,
                      typename RightPhaseSpec,
                      template <typename> class RightCollectionTpl,
                      typename ArgsTmp>
            static ReturnType run(const LeftModelTpl<LeftPhaseSpec, LeftCollectionTpl> &left_model,
                                  const RightModelTpl<RightPhaseSpec, RightCollectionTpl> &right_model,
                                  ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, left_model, right_model);
            }

            // Data + Args
            template <typename LeftPhaseSpec,
                      template <typename> class LeftCollectionTpl,
                      typename RightPhaseSpec,
                      template <typename> class RightCollectionTpl,
                      typename ArgsTmp>
            static ReturnType run(const LeftDataTpl<LeftPhaseSpec, LeftCollectionTpl> &left_data,
                                  const RightDataTpl<RightPhaseSpec, RightCollectionTpl> &right_data,
                                  ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, left_data, right_data);
            }

            // Model only
            template <typename LeftPhaseSpec,
                      template <typename> class LeftCollectionTpl,
                      typename RightPhaseSpec,
                      template <typename> class RightCollectionTpl>
            static ReturnType run(const LeftModelTpl<LeftPhaseSpec, LeftCollectionTpl> &left_model,
                                  const RightModelTpl<RightPhaseSpec, RightCollectionTpl> &right_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, left_model, right_model);
            }

            // Data only
            template <typename LeftPhaseSpec,
                      template <typename> class LeftCollectionTpl,
                      typename RightPhaseSpec,
                      template <typename> class RightCollectionTpl>
            static ReturnType run(const LeftDataTpl<LeftPhaseSpec, LeftCollectionTpl> &left_data,
                                  const RightDataTpl<RightPhaseSpec, RightCollectionTpl> &right_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, left_data, right_data);
            }

            // Base Model + Args
            template <typename LeftModelType, typename RightModelType, typename ArgsTmp>
            static ReturnType run(const LeftModelBaseOf_t<LeftModelType> &left_model,
                                  const RightModelBaseOf_t<RightModelType> &right_model,
                                  ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(left_model.derived(), right_model.derived());
            }

            // Base Data + Args
            template <typename LeftDataType, typename RightDataType, typename ArgsTmp>
            static ReturnType run(const LeftDataBaseOf_t<LeftDataType> &left_data,
                                  const RightDataBaseOf_t<RightDataType> &right_data,
                                  ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(left_data.derived(), right_data.derived());
            }

            // Base Model only
            template <typename LeftModelType, typename RightModelType>
            static ReturnType run(const LeftModelBaseOf_t<LeftModelType> &left_model,
                                  const RightModelBaseOf_t<RightModelType> &right_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(left_model.derived(), right_model.derived());
            }

            // Base Data only
            template <typename LeftDataType, typename RightDataType>
            static ReturnType run(const LeftDataBaseOf_t<LeftDataType> &left_data,
                                  const RightDataBaseOf_t<RightDataType> &right_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(left_data.derived(), right_data.derived());
            }

        private:
            // Internal visitor for Model + Data + Args
            template <typename LeftModel, typename RightModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using LeftData = typename traits<LeftModel>::Data_t;
                using RightData = typename traits<RightModel>::Data_t;

                InternalVisitorModelAndData(LeftData &left_data_, RightData &right_data_, ArgType args_)
                    : left_data(left_data_), right_data(right_data_), args(args_)
                {
                }

                template <typename LeftModelType, typename RightModelType>
                ReturnType operator()(const LeftModelBaseOf_t<LeftModelType> &left_model,
                                      const RightModelBaseOf_t<RightModelType> &right_model) const
                {
                    return bf::invoke(&VisitorDerived::template algo<LeftModelType, RightModelType>,
                                      gf::append(boost::ref(left_model.derived()),
                                                 boost::ref(boost::get<LeftDataOf_t<LeftModelType>>(left_data)),
                                                 boost::ref(boost::get<RightDataOf_t<RightModelType>>(right_data)),
                                                 args));
                }

                LeftData &left_data;
                RightData &right_data;
                ArgType args;
            };

            // Specialization for NoArg
            template <typename LeftModel, typename RightModel>
            struct InternalVisitorModelAndData<LeftModel, RightModel, NoArg> : public boost::static_visitor<ReturnType>
            {
                using LeftData = typename traits<LeftModel>::Data_t;
                using RightData = typename traits<RightModel>::Data_t;

                InternalVisitorModelAndData(LeftData &left_data_, RightData &right_data_)
                    : left_data(left_data_), right_data(right_data_)
                {
                }

                template <typename LeftModelType, typename RightModelType>
                ReturnType operator()(const LeftModelBaseOf_t<LeftModelType> &left_model,
                                      const RightModelBaseOf_t<RightModelType> &right_model) const
                {
                    return bf::invoke(
                        &VisitorDerived::template algo<LeftModelType, RightModelType>,
                        bf::make_vector(boost::ref(left_model.derived()),
                                        boost::ref(boost::get<LeftDataOf_t<LeftModelType>>(left_data)),
                                        boost::ref(boost::get<RightDataOf_t<RightModelType>>(right_data))));
                }

                LeftData &left_data;
                RightData &right_data;
            };

            // Internal visitor for Model/Data only + Args
            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args_) : args(args_) {}

                template <typename LeftModelType, typename RightModelType>
                ReturnType operator()(const LeftModelBaseOf_t<LeftModelType> &left_model,
                                      const RightModelBaseOf_t<RightModelType> &right_model) const
                {
                    return bf::invoke(
                        &VisitorDerived::template algo<LeftModelType, RightModelType>,
                        gf::append(boost::ref(left_model.derived()), boost::ref(right_model.derived()), args));
                }

                template <typename LeftDataType, typename RightDataType>
                ReturnType operator()(const LeftDataBaseOf_t<LeftDataType> &left_data,
                                      const RightDataBaseOf_t<RightDataType> &right_data) const
                {
                    return bf::invoke(
                        &VisitorDerived::template algo<LeftDataType, RightDataType>,
                        gf::append(boost::ref(left_data.derived()), boost::ref(right_data.derived()), args));
                }

                ArgType args;
            };

            // Specialization for NoArg
            template <typename Dummy>
            struct InternalVisitorModel<NoArg, Dummy> : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel() {}

                template <typename LeftModelType, typename RightModelType>
                ReturnType operator()(const LeftModelBaseOf_t<LeftModelType> &left_model,
                                      const RightModelBaseOf_t<RightModelType> &right_model) const
                {
                    return VisitorDerived::template algo<LeftModelType, RightModelType>(left_model.derived(),
                                                                                        right_model.derived());
                }

                template <typename LeftDataType, typename RightDataType>
                ReturnType operator()(const LeftDataBaseOf_t<LeftDataType> &left_data,
                                      const RightDataBaseOf_t<RightDataType> &right_data) const
                {
                    return VisitorDerived::template algo<LeftDataType, RightDataType>(left_data.derived(),
                                                                                      right_data.derived());
                }
            };

        }; // struct BinaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_common_visitors_binary_visitor_hpp__
