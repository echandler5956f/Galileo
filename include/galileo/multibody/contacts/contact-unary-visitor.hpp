#ifndef __galileo_multibody_contacts_contact_unary_visitor_hpp__
#define __galileo_multibody_contacts_contact_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/common/meta/fusion.hpp"
#include "galileo/multibody/contacts/contact-base.hpp"

namespace galileo
{
    namespace fusion
    {

        // Base structure for Unary visitation of a ContactModel.
        // This structure provides runners to call the right visitor according to the number of
        // arguments.
        template <typename ContactVisitorDerived, typename ReturnType = void>
        struct ContactUnaryVisitorBase
        {
            template <
                typename PhaseSpec,
                template <typename PS> class ContactCollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
                ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ContactModelTpl<PhaseSpec, ContactCollectionTpl>, ArgsTmp>
                    visitor(contact_data, args);
                return boost::apply_visitor(visitor, contact_model);
            }

            template <typename PhaseSpec, template <typename PS> class ContactCollectionTpl>
            static ReturnType run(
                const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
                ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
            {
                InternalVisitorModelAndData<ContactModelTpl<PhaseSpec, ContactCollectionTpl>, NoArg>
                    visitor(contact_data);
                return boost::apply_visitor(visitor, contact_model);
            }

            template <typename ContactModelType, typename ArgsTmp>
            static ReturnType run(
                const ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS> &contact_model,
                typename ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS>::Data_t &contact_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ContactModelType, ArgsTmp> visitor(contact_data, args);
                return visitor(contact_model.derived());
            }

            template <typename ContactModelType>
            static ReturnType run(
                const ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS> &contact_model,
                typename ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS>::Data_t &contact_data)
            {
                InternalVisitorModelAndData<ContactModelType, NoArg> visitor(contact_data);
                return visitor(contact_model.derived());
            }

            template <
                typename PhaseSpec,
                template <typename PS> class ContactCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, contact_model);
            }

            template <
                typename PhaseSpec,
                template <typename PS> class ContactCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, contact_data);
            }

            template <typename PhaseSpec, template <typename PS> class ContactCollectionTpl>
            static ReturnType run(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, contact_model);
            }

            template <typename PhaseSpec, template <typename PS> class ContactCollectionTpl>
            static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, contact_data);
            }

            template <typename ContactModelType, typename ArgsTmp>
            static ReturnType run(const ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS> &contact_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(contact_model.derived());
            }

            template <typename ContactDataType, typename ArgsTmp>
            static ReturnType run(const ContactDataBase<ContactDataType, typename traits<ContactDataType>::PS> &contact_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(contact_data.derived());
            }

            template <typename ContactModelType>
            static ReturnType run(const ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS> &contact_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(contact_model.derived());
            }

            template <typename ContactDataType>
            static ReturnType run(const ContactDataBase<ContactDataType, typename traits<ContactDataType>::PS> &contact_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(contact_data.derived());
            }

        private:
            template <typename ContactModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using ContactData = typename traits<ContactModel>::Data_t;

                InternalVisitorModelAndData(ContactData &contact_data, ArgType args)
                    : contact_data(contact_data), args(args)
                {
                }

                template <typename ContactModelType>
                ReturnType operator()(const ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS> &contact_model) const
                {
                    return bf::invoke(
                        &ContactVisitorDerived::template algo<ContactModelType>,
                        gf::append(
                            boost::ref(contact_model.derived()),
                            boost::ref(
                                boost::get<typename ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS>::Data_t>(contact_data)),
                            args));
                }

                ReturnType operator()(const ContactModelVoid)
                {
                    return;
                }

                ContactData &contact_data;
                ArgType args;
            };

            template <typename ContactModel>
            struct InternalVisitorModelAndData<ContactModel, NoArg>
                : public boost::static_visitor<ReturnType>
            {
                using ContactData = typename traits<ContactModel>::Data_t;

                InternalVisitorModelAndData(ContactData &contact_data)
                    : contact_data(contact_data)
                {
                }

                template <typename ContactModelType>
                ReturnType operator()(const ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS> &contact_model) const
                {
                    return bf::invoke(
                        &ContactVisitorDerived::template algo<ContactModelType>,
                        bf::make_vector(
                            boost::ref(contact_model.derived()),
                            boost::ref(
                                boost::get<typename ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS>::Data_t>(contact_data))));
                }

                ContactData &contact_data;
            };

            template <typename ArgType, typename Dummy = void>
            struct InternalVisitorModel : public boost::static_visitor<ReturnType>
            {
                InternalVisitorModel(ArgType args)
                    : args(args)
                {
                }

                template <typename ContactModelType>
                ReturnType operator()(const ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS> &contact_model) const
                {
                    return bf::invoke(
                        &ContactVisitorDerived::template algo<ContactModelType>,
                        gf::append(boost::ref(contact_model.derived()), args));
                }

                template <typename ContactDataType>
                ReturnType operator()(const ContactDataBase<ContactDataType, typename traits<ContactDataType>::PS> &contact_data) const
                {
                    return bf::invoke(
                        &ContactVisitorDerived::template algo<ContactDataType>,
                        gf::append(boost::ref(contact_data.derived()), args));
                }

                ReturnType operator()(const ContactModelVoid)
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

                template <typename ContactModelType>
                ReturnType operator()(const ContactModelBase<ContactModelType, typename traits<ContactModelType>::PS> &contact_model) const
                {
                    return ContactVisitorDerived::template algo<ContactModelType>(contact_model.derived());
                }

                template <typename ContactDataType>
                ReturnType operator()(const ContactDataBase<ContactDataType, typename traits<ContactDataType>::PS> &contact_data) const
                {
                    return ContactVisitorDerived::template algo<ContactDataType>(contact_data.derived());
                }
            };

        }; // struct ContactUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_unary_visitor_hpp__
