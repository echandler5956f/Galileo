#ifndef __galileo_multibody_contacts_contact_unary_visitor_hpp__
#define __galileo_multibody_contacts_contact_unary_visitor_hpp__

#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/get.hpp>

#include "galileo/utils/fusion.hpp"
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
                const galileo::ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
                galileo::ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<galileo::ContactModelTpl<PhaseSpec, ContactCollectionTpl>, ArgsTmp>
                    visitor(contact_data, args);
                return boost::apply_visitor(visitor, contact_model);
            }

            template <typename PhaseSpec, template <typename PS> class ContactCollectionTpl>
            static ReturnType run(
                const galileo::ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
                galileo::ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
            {
                InternalVisitorModelAndData<galileo::ContactModelTpl<PhaseSpec, ContactCollectionTpl>, NoArg>
                    visitor(contact_data);
                return boost::apply_visitor(visitor, contact_model);
            }

            template <typename ContactModelDerived, typename ArgsTmp>
            static ReturnType run(
                const galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS> &contact_model,
                typename galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS>::ContactDataDerived &contact_data,
                ArgsTmp args)
            {
                InternalVisitorModelAndData<ContactModelDerived, ArgsTmp> visitor(contact_data, args);
                return visitor(contact_model.derived());
            }

            template <typename ContactModelDerived>
            static ReturnType run(
                const galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS> &contact_model,
                typename galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS>::ContactDataDerived &contact_data)
            {
                InternalVisitorModelAndData<ContactModelDerived, NoArg> visitor(contact_data);
                return visitor(contact_model.derived());
            }

            template <
                typename PhaseSpec,
                template <typename PS> class ContactCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const galileo::ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, contact_model);
            }

            template <
                typename PhaseSpec,
                template <typename, int> class ContactCollectionTpl,
                typename ArgsTmp>
            static ReturnType
            run(const galileo::ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return boost::apply_visitor(visitor, contact_data);
            }

            template <typename PhaseSpec, template <typename PS> class ContactCollectionTpl>
            static ReturnType run(const galileo::ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, contact_model);
            }

            template <typename PhaseSpec, template <typename PS> class ContactCollectionTpl>
            static ReturnType run(const galileo::ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return boost::apply_visitor(visitor, contact_data);
            }

            template <typename ContactModelDerived, typename ArgsTmp>
            static ReturnType run(const galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS> &contact_model, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(contact_model.derived());
            }

            template <typename ContactDataDerived, typename ArgsTmp>
            static ReturnType run(const galileo::ContactDataBase<ContactDataDerived, typename traits<ContactDataDerived>::PS> &contact_data, ArgsTmp args)
            {
                InternalVisitorModel<ArgsTmp> visitor(args);
                return visitor(contact_data.derived());
            }

            template <typename ContactModelDerived>
            static ReturnType run(const galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS> &contact_model)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(contact_model.derived());
            }

            template <typename ContactDataDerived>
            static ReturnType run(const galileo::ContactDataBase<ContactDataDerived, typename traits<ContactDataDerived>::PS> &contact_data)
            {
                InternalVisitorModel<NoArg> visitor;
                return visitor(contact_data.derived());
            }

        private:
            template <typename ContactModel, typename ArgType>
            struct InternalVisitorModelAndData : public boost::static_visitor<ReturnType>
            {
                using ContactData = typename traits<ContactModel>::ContactDataDerived;

                InternalVisitorModelAndData(ContactData &contact_data, ArgType args)
                    : contact_data(contact_data), args(args)
                {
                }

                template <typename ContactModelDerived>
                ReturnType operator()(const galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS> &contact_model) const
                {
                    return bf::invoke(
                        &ContactVisitorDerived::template algo<ContactModelDerived>,
                        bf::append(
                            boost::ref(contact_model.derived()),
                            boost::ref(
                                boost::get<typename galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS>::ContactDataDerived>(contact_data)),
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
                using ContactData = typename traits<ContactModel>::ContactDataDerived;

                InternalVisitorModelAndData(ContactData &contact_data)
                    : contact_data(contact_data)
                {
                }

                template <typename ContactModelDerived>
                ReturnType operator()(const galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS> &contact_model) const
                {
                    return bf::invoke(
                        &ContactVisitorDerived::template algo<ContactModelDerived>,
                        bf::make_vector(
                            boost::ref(contact_model.derived()),
                            boost::ref(
                                boost::get<typename galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS>::ContactDataDerived>(contact_data))));
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

                template <typename ContactModelDerived>
                ReturnType operator()(const galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS> &contact_model) const
                {
                    return bf::invoke(
                        &ContactVisitorDerived::template algo<ContactModelDerived>,
                        bf::append(boost::ref(contact_model.derived()), args));
                }

                template <typename ContactDataDerived>
                ReturnType operator()(const galileo::ContactDataBase<ContactDataDerived, typename traits<ContactDataDerived>::PS> &contact_data) const
                {
                    return bf::invoke(
                        &ContactVisitorDerived::template algo<ContactDataDerived>,
                        bf::append(boost::ref(contact_data.derived()), args));
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

                template <typename ContactModelDerived>
                ReturnType operator()(const galileo::ContactModelBase<ContactModelDerived, typename traits<ContactModelDerived>::PS> &contact_model) const
                {
                    return ContactVisitorDerived::template algo<ContactModelDerived>(contact_model.derived());
                }

                template <typename ContactDataDerived>
                ReturnType operator()(const galileo::ContactDataBase<ContactDataDerived, typename traits<ContactDataDerived>::PS> &contact_data) const
                {
                    return ContactVisitorDerived::template algo<ContactDataDerived>(contact_data.derived());
                }
            };

        }; // struct ContactUnaryVisitorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_unary_visitor_hpp__