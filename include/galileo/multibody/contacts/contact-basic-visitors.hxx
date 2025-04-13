#ifndef __galileo_multibody_contacts_contact_basic_visitors_hxx__
#define __galileo_multibody_contacts_contact_basic_visitors_hxx__

#include <vector>

#include <boost/fusion/container/generation/make_vector.hpp>
#include "galileo/multibody/contacts/contact-unary-visitor.hpp"

#include "galileo/multibody/contacts/contact-basic-visitors.hpp"

#include "galileo/utils/aligned-vector.hpp"

namespace galileo
{

    // Contact model visitors

    template <typename PhaseSpec,
              template <typename> class ContactCollectionTpl,
              typename DataCollector>
    struct ContactCreateDataVisitor
        : fusion::ContactUnaryVisitorBase<ContactCreateDataVisitor<PhaseSpec, ContactCollectionTpl, DataCollector>>
    {
        using ArgsType = boost::fusion::vector<DataCollector *const>;
        using ContactCollection_t = ContactCollectionTpl<PhaseSpec>;
        using ContactModelVariant_t = ContactCollection_t::ModelVariant_t;
        using ContactDataVariant_t = ContactDataTpl<PhaseSpec, ContactCollectionTpl>;

        template <typename ContactModelDerived>
        static ContactDataVariant_t algo(
            const ContactModelBase<ContactModelDerived, PhaseSpec> &contact_model,
            DataCollector *const collector)
        {
            return ContactDataVariant_t(contact_model.createData(collector));
        }
    };

    template <typename PhaseSpec,
              template <typename> class ContactCollectionTpl,
              typename DataCollector>
    inline ContactDataTpl<PhaseSpec, ContactCollectionTpl> contact_create_data(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        DataCollector *const collector)
    {
        typedef ContactCreateDataVisitor<PhaseSpec, ContactCollectionTpl, DataCollector> Algo;

        return Algo::run(contact_model, typename Algo::ArgsType(collector));
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename StateVectorType>
    struct ContactCalcZerothOrderVisitor
        : fusion::ContactUnaryVisitorBase<ContactCalcZerothOrderVisitor<PhaseSpec, ContactCollectionTpl, StateVectorType>>
    {
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>>;

        template <typename ContactModel>
        static void algo(
            const ContactModelBase<ContactModel, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModel::ContactDataDerived, PhaseSpec> &contact_data,
            const Eigen::MatrixBase<StateVectorType> &x)
        {
            contact_model.calc(contact_data, x.derived());
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
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>>;

        template <typename ContactModel>
        static void algo(
            const ContactModelBase<ContactModel, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModel::ContactDataDerived, PhaseSpec> &contact_data,
            const Eigen::MatrixBase<StateVectorType> &x)
        {
            contact_model.calcDiff(contact_data, x.derived());
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
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<ForceVectorType>>;

        template <typename ContactModel>
        static void algo(
            const ContactModelBase<ContactModel, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModel::ContactDataDerived, PhaseSpec> &contact_data,
            const Eigen::MatrixBase<ForceVectorType> &force)
        {
            contact_model.updateForce(contact_data, force.derived());
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
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<MatrixNcNdxType>, Eigen::MatrixBase<MatrixNcNuType>>;

        template <typename ContactModel>
        static void algo(
            const ContactModelBase<ContactModel, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModel::ContactDataDerived, PhaseSpec> &contact_data,
            const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
            const Eigen::MatrixBase<MatrixNcNuType> &df_du)
        {
            contact_model.updateForceDiff(contact_data, df_dx.derived(), df_du.derived());
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

        template <typename ContactModel>
        static void algo(
            const ContactModelBase<ContactModel, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModel::ContactDataDerived, PhaseSpec> &contact_data)
        {
            contact_model.setZeroForce(contact_data);
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

        template <typename ContactModel>
        static void algo(
            const ContactModelBase<ContactModel, PhaseSpec> &contact_model,
            ContactDataBase<typename ContactModel::ContactDataDerived, PhaseSpec> &contact_data)
        {
            contact_model.setZeroForceDiff(contact_data);
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
              template <typename PS> class ContactCollectionTpl>
    struct ContactNcVisitor : boost::static_visitor<int>
    {
        template <typename ContactModelDerived>
        int operator()(const ContactModelBase<ContactModelDerived, PhaseSpec> &contact_model) const
        {
            return contact_model.nc();
        }

        static int run(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
        {
            return boost::apply_visitor(ContactNcVisitor<PhaseSpec, ContactCollectionTpl>(), contact_model);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline int contact_nc(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
    {
        return ContactNcVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_model);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactIdVisitor : boost::static_visitor<int>
    {
        using ReturnType = typename traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>::Index_t;

        template <typename ContactModelDerived>
        ReturnType operator()(const ContactModelBase<ContactModelDerived, PhaseSpec> &contact_model) const
        {
            return contact_model.id();
        }

        static ReturnType run(const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
        {
            return boost::apply_visitor(ContactIdVisitor<PhaseSpec, ContactCollectionTpl>(), contact_model);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>::Index_t contact_id(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model)
    {
        return ContactIdVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_model);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactSetIdVisitor : boost::static_visitor<void>
    {
        using ArgsType = boost::fusion::vector<typename traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>::Index_t>;

        template <typename ContactModel>
        static void algo(
            const ContactModelBase<ContactModel, PhaseSpec> &contact_model,
            const typename traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>::Index_t &id)
        {
            contact_model.set_id(id);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline void contact_set_id(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        const typename traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>::Index_t &id)
    {
        typedef ContactSetIdVisitor<PhaseSpec, ContactCollectionTpl> Algo;

        Algo::run(contact_model, typename Algo::ArgsType(id));
    }

    // Contact data visitors

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactRobotDataPointerVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::RobotDataPointer_t>
    {

        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::RobotDataPointer_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
        {
            return contact_data.robot_data_pointer();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactRobotDataPointerVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::RobotDataPointer_t contact_robot_data_pointer(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactRobotDataPointerVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactFrameVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Index_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Index_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Index_t contact_frame(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactFrameVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactTypeVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::ReferenceFrame_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::ReferenceFrame_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
        {
            return contact_data.type();
        }

        static ReturnType run(const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
        {
            return boost::apply_visitor(ContactTypeVisitor<PhaseSpec, ContactCollectionTpl>(), contact_data);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::ReferenceFrame_t contact_type(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactTypeVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactJmFVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::SE3_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::SE3_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::SE3_t contact_jMf(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactJmFVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactJcVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNv_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNv_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNv_t contact_Jc(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactJcVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactFVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Force_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Force_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Force_t contact_f(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactFVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactFextVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Force_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Force_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Force_t contact_fext(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactFextVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactdFdXVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNdx_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNdx_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNdx_t contact_df_dx(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactdFdXVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactdFdUVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNu_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNu_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNu_t contact_df_du(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactdFdUVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactFXjVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::ActionMatrix_t>
    {

        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::ActionMatrix_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::ActionMatrix_t contact_fXj(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactFXjVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactA0Visitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::VectorNc_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::VectorNc_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::VectorNc_t contact_a0(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactA0Visitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactDA0dXVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNdx_t>
    {

        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNdx_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNdx_t contact_da0_dx(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactDA0dXVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactDtauDqVisitor : boost::static_visitor<typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNv_t>
    {
        using ReturnType = typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNv_t;

        template <typename ContactDataDerived>
        ReturnType operator()(const ContactDataBase<ContactDataDerived, PhaseSpec> &contact_data) const
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
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNv_t contact_dtau_dq(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data)
    {
        return ContactDtauDqVisitor<PhaseSpec, ContactCollectionTpl>::run(contact_data);
    }

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_basic_visitors_hxx__