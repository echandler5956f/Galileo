#ifndef __galileo_multibody_contacts_contact_generic_hpp__
#define __galileo_multibody_contacts_contact_generic_hpp__

#include "galileo/multibody/contacts/fwd.hpp"
#include "galileo/multibody/contacts/contact-base.hpp"
#include "galileo/multibody/contacts/contact-collection.hpp"
#include "galileo/multibody/contacts/contact-basic-visitors.hxx"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactTpl;

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<ContactTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ContactTpl<PS, ContactCollectionTpl>;
        using Collection_t = ContactCollectionTpl<PS>;
        using Model_t = ContactModelTpl<PS, ContactCollectionTpl>;
        using Data_t = ContactDataTpl<PS, ContactCollectionTpl>;

        static constexpr int NC = Eigen::Dynamic;

        // Traits required by ForceDataBase
        using RobotData_t = typename PS::RobotData_t;
        using Index_t = typename PS::Index_t;
        using ReferenceFrame_t = typename PS::ReferenceFrame_t;
        using SE3_t = typename PS::SE3_t;
        using MatrixNcNv_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NV, PS::Options, 6, PS::NV>;
        using Force_t = typename PS::Force_t;
        using MatrixNcNdx_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NDX, PS::Options, 6, PS::NDX>;
        using MatrixNcNu_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NU, PS::Options, 6, PS::NU>;

        // Traits required by ContactDataBase
        using ActionMatrix_t = typename PS::ActionMatrix_t;
        using VectorNc_t = Eigen::Matrix<typename PS::VarScalar, NC, 1, PS::Options, 6, 1>;
        using MatrixNv_t = Eigen::Matrix<typename PS::VarScalar, PS::NV, PS::NV, PS::Options>;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ContactTpl<PS, ContactCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ContactTpl<PS, ContactCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactDataTpl : public ContactDataBase<ContactDataTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                            ContactCollectionTpl<PhaseSpec>::DataVariant_t
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ContactTpl<PS, ContactCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_FORCE_DATA_TYPEDEF(Meta_t);
        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);

        using DataVariant_t = typename Collection_t::DataVariant_t;

        DataVariant_t &toVariant()
        {
            return *static_cast<DataVariant_t *>(this);
        }
        const DataVariant_t &toVariant() const
        {
            return *static_cast<const DataVariant_t *>(this);
        }

        RobotData_t *robot_data_pointer() const
        {
            return galileo::contact_robot_data_pointer(*this);
        }

        Index_t frame() const
        {
            return galileo::contact_frame(*this);
        }

        ReferenceFrame_t type() const
        {
            return galileo::contact_type(*this);
        }

        SE3_t jMf() const
        {
            return galileo::contact_jMf(*this);
        }

        MatrixNcNv_t Jc() const
        {
            return galileo::contact_Jc(*this);
        }

        Force_t f() const
        {
            return galileo::contact_f(*this);
        }

        Force_t fext() const
        {
            return galileo::contact_fext(*this);
        }

        MatrixNcNdx_t df_dx() const
        {
            return galileo::contact_df_dx(*this);
        }

        MatrixNcNu_t df_du() const
        {
            return galileo::contact_df_du(*this);
        }

        ActionMatrix_t fXj() const
        {
            return galileo::contact_fXj(*this);
        }

        VectorNc_t a0() const
        {
            return galileo::contact_a0(*this);
        }

        MatrixNcNdx_t da0_dx() const
        {
            return galileo::contact_da0_dx(*this);
        }

        MatrixNv_t dtau_dq() const
        {
            return galileo::contact_dtau_dq(*this);
        }

        ContactDataTpl()
            : DataVariant_t()
        {
        }

        ContactDataTpl(const DataVariant_t &data_variant)
            : DataVariant_t(data_variant)
        {
        }

        template <typename DataDerived>
        ContactDataTpl(const ContactDataBase<DataDerived, PhaseSpec> &data)
            : Collection_t::DataVariant_t((DataVariant_t)data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename DataVariant_t::types, DataDerived>));
        }

        GENERIC_ACCESSOR(RobotData_t *, robot_data_pointer);
        GENERIC_ACCESSOR(Index_t, frame);
        GENERIC_ACCESSOR(ReferenceFrame_t, type);
        GENERIC_ACCESSOR(SE3_t, jMf);
        GENERIC_ACCESSOR(MatrixNcNv_t, Jc);
        GENERIC_ACCESSOR(Force_t, f);
        GENERIC_ACCESSOR(Force_t, fext);
        GENERIC_ACCESSOR(MatrixNcNdx_t, df_dx);
        GENERIC_ACCESSOR(MatrixNcNu_t, df_du);

        GENERIC_ACCESSOR(ActionMatrix_t, fXj);
        GENERIC_ACCESSOR(VectorNc_t, a0);
        GENERIC_ACCESSOR(MatrixNcNdx_t, da0_dx);
        GENERIC_ACCESSOR(MatrixNv_t, dtau_dq);
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactModelTpl : public ContactModelBase<ContactModelTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                             ContactCollectionTpl<PhaseSpec>::ModelVariant_t
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ContactTpl<PS, ContactCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using ModelVariant_t = typename Collection_t::ModelVariant_t;

        using Index_t = typename traits<Meta_t>::Index_t;

        ContactModelTpl()
            : ModelVariant_t()
        {
        }

        ContactModelTpl(const ModelVariant_t &model_variant)
            : ModelVariant_t(model_variant)
        {
        }

        template <typename ModelDerived>
        ContactModelTpl(const ContactModelBase<ModelDerived, PhaseSpec> &model)
            : Collection_t::ModelVariant_t((ModelVariant_t)model.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename ModelVariant_t::types, ModelDerived>));
        }

        ModelVariant_t &toVariant()
        {
            return *static_cast<ModelVariant_t *>(this);
        }

        const ModelVariant_t &toVariant() const
        {
            return *static_cast<const ModelVariant_t *>(this);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector)
        {
            return galileo::contact_create_data(*this, collector);
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x)
        {
            galileo::contact_calc_zeroth_order(*this, data, x);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
        {
            galileo::contact_calc_first_order(*this, data, x);
        }

        template <typename ForceVectorType>
        void updateForce(Data_t &data,
                         const Eigen::MatrixBase<ForceVectorType> &force)
        {
            galileo::contact_update_force(*this, data, force.derived());
        }

        template <typename MatrixNcNdxType, typename MatrixNcNuType>
        void updateForceDiff(Data_t &data,
                             const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
                             const Eigen::MatrixBase<MatrixNcNuType> &df_du) const
        {
            galileo::contact_update_force_diff(*this, data, df_dx.derived(), df_du.derived());
        }

        void setZeroForce(Data_t &data) const
        {
            galileo::contact_set_zero_force(*this, data);
        }

        void setZeroForceDiff(Data_t &data) const
        {
            galileo::contact_set_zero_force_diff(*this, data);
        }

        int nc_impl() const
        {
            return galileo::contact_nc(*this);
        }

        Index_t id_impl() const
        {
            return galileo::contact_id(*this);
        }

        void set_id_impl(const Index_t &id)
        {
            galileo::contact_set_id(*this, id);
        }

    }; // struct ContactModelTpl

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_generic_hpp__
