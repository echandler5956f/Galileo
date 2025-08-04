#ifndef __galileo_multibody_contacts_contact_visitors_hxx__
#define __galileo_multibody_contacts_contact_visitors_hxx__

#include "galileo/multibody/contacts/contact-visitor-base.hpp"
#include "galileo/multibody/contacts/contact-visitors.hpp"

namespace galileo
{

    // Contact model visitors

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename StateVectorType>
    struct ContactCalcZerothOrderVisitor
        : fusion::ContactUnaryVisitorBase<ContactCalcZerothOrderVisitor<PhaseSpec, ContactCollectionTpl, StateVectorType>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &>;

        template <typename ContactModelType>
        static void algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModelType::Data_t, PhaseSpec> &contact_data,
            const Eigen::MatrixBase<StateVectorType> &x)
        {
            contact_model.calc(contact_data.derived(), x.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename StateVectorType>
    inline void contact_calc_zeroth_order(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
        const Eigen::MatrixBase<StateVectorType> &x)
    {
        typedef ContactCalcZerothOrderVisitor<PhaseSpec, ContactCollectionTpl, StateVectorType> Algo;

        Algo::run(contact_model, contact_data, typename Algo::ArgsType(x.derived()));
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename StateVectorType>
    struct ContactCalcFirstOrderVisitor
        : fusion::ContactUnaryVisitorBase<ContactCalcFirstOrderVisitor<PhaseSpec, ContactCollectionTpl, StateVectorType>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &>;

        template <typename ContactModelType>
        static void algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModelType::Data_t, PhaseSpec> &contact_data,
            const Eigen::MatrixBase<StateVectorType> &x)
        {
            contact_model.calcDiff(contact_data.derived(), x.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename StateVectorType>
    inline void contact_calc_first_order(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
        const Eigen::MatrixBase<StateVectorType> &x)
    {
        typedef ContactCalcFirstOrderVisitor<PhaseSpec, ContactCollectionTpl, StateVectorType> Algo;

        Algo::run(contact_model, contact_data, typename Algo::ArgsType(x.derived()));
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename ForceVectorType>
    struct ContactUpdateForceVisitor
        : fusion::ContactUnaryVisitorBase<ContactUpdateForceVisitor<PhaseSpec, ContactCollectionTpl, ForceVectorType>>
    {
        using ArgsType = boost::fusion::vector<const ForceVectorType &>;

        template <typename ContactModelType>
        static void algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModelType::Data_t, PhaseSpec> &contact_data,
            const Eigen::MatrixBase<ForceVectorType> &force)
        {
            contact_model.updateForce(contact_data.derived(), force.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename ForceVectorType>
    inline void contact_update_force(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
        const Eigen::MatrixBase<ForceVectorType> &force)
    {
        typedef ContactUpdateForceVisitor<PhaseSpec, ContactCollectionTpl, ForceVectorType> Algo;

        Algo::run(contact_model, contact_data, typename Algo::ArgsType(force.derived()));
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename MatrixNcNdxType,
              typename MatrixNcNuType>
    struct ContactUpdateForceDiffVisitor
        : fusion::ContactUnaryVisitorBase<ContactUpdateForceDiffVisitor<PhaseSpec, ContactCollectionTpl, MatrixNcNdxType, MatrixNcNuType>>
    {
        using ArgsType = boost::fusion::vector<const MatrixNcNdxType &, const MatrixNcNuType &>;

        template <typename ContactModelType>
        static void algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModelType::Data_t, PhaseSpec> &contact_data,
            const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
            const Eigen::MatrixBase<MatrixNcNuType> &df_du)
        {
            contact_model.updateForceDiff(contact_data.derived(), df_dx.derived(), df_du.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename MatrixNcNdxType,
              typename MatrixNcNuType>
    inline void contact_update_force_diff(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
        const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
        const Eigen::MatrixBase<MatrixNcNuType> &df_du)
    {
        typedef ContactUpdateForceDiffVisitor<PhaseSpec, ContactCollectionTpl, MatrixNcNdxType, MatrixNcNuType> Algo;

        Algo::run(contact_model, contact_data, typename Algo::ArgsType(df_dx.derived(), df_du.derived()));
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactSetZeroForceVisitor
        : fusion::ContactUnaryVisitorBase<ContactSetZeroForceVisitor<PhaseSpec, ContactCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<>;

        template <typename ContactModelType>
        static void algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModelType::Data_t, PhaseSpec> &contact_data)
        {
            contact_model.setZeroForce(contact_data.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline void contact_set_zero_force(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        typedef ContactSetZeroForceVisitor<PhaseSpec, ContactCollectionTpl> Algo;

        Algo::run(contact_model, contact_data, typename Algo::ArgsType());
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactSetZeroForceDiffVisitor
        : fusion::ContactUnaryVisitorBase<ContactSetZeroForceDiffVisitor<PhaseSpec, ContactCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<>;

        template <typename ContactModelType>
        static void algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModelType::Data_t, PhaseSpec> &contact_data)
        {
            contact_model.setZeroForceDiff(contact_data.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline void contact_set_zero_force_diff(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        typedef ContactSetZeroForceDiffVisitor<PhaseSpec, ContactCollectionTpl> Algo;

        Algo::run(contact_model, contact_data, typename Algo::ArgsType());
    }

    template <typename PhaseSpec,
              template <typename> class ContactCollectionTpl>
    struct ContactCreateDataVisitor
        : fusion::ContactUnaryVisitorBase<ContactCreateDataVisitor<PhaseSpec, ContactCollectionTpl>,
                                          ContactDataTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<typename PhaseSpec::RobotData_t *const>;
        using ContactCollection_t = ContactCollectionTpl<PhaseSpec>;
        using ContactModelVariant_t = ContactCollection_t::ContactModelVariant_t;
        using ContactDataVariant_t = ContactDataTpl<PhaseSpec, ContactCollectionTpl>;

        template <typename ContactModelType>
        static ContactDataVariant_t algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            typename PhaseSpec::RobotData_t *const robot)
        {
            return ContactDataVariant_t(contact_model.createData(robot));
        }
    };

    template <typename PhaseSpec,
              template <typename> class ContactCollectionTpl>
    inline ContactDataTpl<PhaseSpec, ContactCollectionTpl> contact_create_data(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        typename PhaseSpec::RobotData_t *const robot)
    {
        typedef ContactCreateDataVisitor<PhaseSpec, ContactCollectionTpl> Algo;
        return Algo::run(contact_model, typename Algo::ArgsType(robot));
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactGetIdVisitor : boost::static_visitor<typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t>
    {
        using ReturnType = typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t;

        template <typename ContactModelType>
        ReturnType operator()(const ContactModelBase<ContactModelType, PhaseSpec> &contact_model) const
        {
            return contact_model.get_id();
        }

        static ReturnType run(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
        {
            return boost::apply_visitor(ContactGetIdVisitor<PhaseSpec, ContactCollectionTpl>(), contact_model);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t contact_get_id(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
    {
        return ContactGetIdVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_model);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactSetIdVisitor : boost::static_visitor<void>
    {
        using ArgsType = boost::fusion::vector<const typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t &>;

        template <typename ContactModelType>
        static void algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            const typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t &id)
        {
            contact_model.set_id(id);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline void contact_set_id(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        const typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t &id)
    {
        typedef ContactSetIdVisitor<PhaseSpec, ContactCollectionTpl> Algo;

        Algo::run(contact_model, typename Algo::ArgsType(id));
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactGetTypeVisitor : boost::static_visitor<typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t>
    {
        using ReturnType = typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t;

        template <typename ContactModelType>
        ReturnType operator()(const ContactModelBase<ContactModelType, PhaseSpec> &contact_model) const
        {
            return contact_model.get_type();
        }

        static ReturnType run(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
        {
            return boost::apply_visitor(ContactGetTypeVisitor<PhaseSpec, ContactCollectionTpl>(), contact_model);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t contact_get_type(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
    {
        return ContactGetTypeVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_model);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactSetTypeVisitor : boost::static_visitor<void>
    {
        using ArgsType = boost::fusion::vector<const typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t &>;

        template <typename ContactModelType>
        static void algo(
            const ContactModelBase<ContactModelType, PhaseSpec> &contact_model,
            const typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t &type)
        {
            contact_model.set_type(type);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline void contact_set_type(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        const typename ContactModelTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t &type)
    {
        typedef ContactSetTypeVisitor<PhaseSpec, ContactCollectionTpl> Algo;

        Algo::run(contact_model, typename Algo::ArgsType(type));
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactGetNcVisitor : boost::static_visitor<int>
    {
        template <typename ContactModelType>
        int operator()(const ContactModelBase<ContactModelType, PhaseSpec> &contact_model) const
        {
            return contact_model.get_nc();
        }

        static int run(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
        {
            return boost::apply_visitor(ContactGetNcVisitor<PhaseSpec, ContactCollectionTpl>(), contact_model);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline int contact_get_nc(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
    {
        return ContactGetNcVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_model);
    }

    // Contact data visitors

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactRobotDataVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::RobotData_t>
    {

        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::RobotDataPointer_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.robot();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactRobotDataVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::RobotData_t *contact_robot_data(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactRobotDataVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactFrameVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.frame();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactFrameVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::FrameIndex_t contact_frame(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactFrameVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactTypeDataVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.type();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactTypeDataVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::ReferenceFrame_t contact_type_data(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactTypeDataVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactJmFVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::SE3_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::SE3_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.jMf();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactJmFVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::SE3_t contact_jMf(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactJmFVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactJcVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNv_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNv_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.Jc();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactJcVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNv_t contact_Jc(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactJcVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactFVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::Force_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::Force_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.f();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactFVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::Force_t contact_f(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactFVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactFextVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::Force_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::Force_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.fext();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactFextVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::Force_t contact_fext(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactFextVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactdFdXVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNdx_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNdx_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.df_dx();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactdFdXVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNdx_t contact_df_dx(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactdFdXVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactdFdUVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNu_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNu_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.df_du();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactdFdUVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNu_t contact_df_du(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactdFdUVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactFXjVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::ActionMatrix_t>
    {

        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::ActionMatrix_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.fXj();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactFXjVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::ActionMatrix_t contact_fXj(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactFXjVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactA0Visitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::VectorNc_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::VectorNc_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.a0();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactA0Visitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::VectorNc_t contact_a0(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactA0Visitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactDA0dXVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNdx_t>
    {

        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNdx_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.da0_dx();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactDA0dXVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNcNdx_t contact_da0_dx(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactDA0dXVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactDtauDqVisitor : boost::static_visitor<typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNv_t>
    {
        using ReturnType = typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNv_t;

        template <typename ContactDataType>
        ReturnType operator()(const ContactDataBase<ContactDataType, PhaseSpec> &contact_data) const
        {
            return contact_data.dtau_dq();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactDtauDqVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename ContactDataTpl<PhaseSpec, ContactCollectionTpl>::MatrixNv_t contact_dtau_dq(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactDtauDqVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_visitors_hxx__
