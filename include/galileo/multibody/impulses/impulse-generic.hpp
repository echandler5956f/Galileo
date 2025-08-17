#ifndef __galileo_multibody_impulses_impulse_generic_hpp__
#define __galileo_multibody_impulses_impulse_generic_hpp__

#include "galileo/multibody/impulses/impulse-base.hpp"
#include "galileo/multibody/impulses/impulse-collection.hpp"
#include "galileo/multibody/impulses/impulse-visitors.hxx"

namespace galileo
{

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseTpl;

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct traits<ImpulseTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = ImpulseTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = ImpulseCollectionTpl<PS>;
        using Model_t = ImpulseModelTpl<PS, ImpulseCollectionTpl>;
        using Data_t = ImpulseDataTpl<PS, ImpulseCollectionTpl>;

        using DimNC_t = DimensionTpl<>;
        static constexpr int NC = DimNC_t::Value;

        using DimNU_t = traits<typename PS::NodeMeta_t>::DimNU_t;

        // Traits required by ForceDataBase
        using MatrixNcNv_t =
            Eigen::GMatrix<typename PS::VarScalar, NC, PS::DimNV_t::Value, PS::Options, 6, PS::DimNV_t::Value>;
        using MatrixNcNdx_t =
            Eigen::GMatrix<typename PS::VarScalar, NC, PS::DimNDX_t::Value, PS::Options, 6, PS::DimNDX_t::Value>;
        using MatrixNcNu_t = Eigen::GMatrix<typename PS::VarScalar, NC, DimNU_t::Value, PS::Options, 6, DimNU_t::Value>;

        // Traits required by ImpulseDataBase
        using VectorNc_t = Eigen::GMatrix<typename PS::VarScalar, NC, 1, PS::Options, 6, 1>;
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct traits<ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = ImpulseTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_IMPULSE_DATA_TYPEDEF(Meta_t);
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct traits<ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = ImpulseTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseDataTpl : public ImpulseDataBase<ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>, PhaseSpec>,
                            ImpulseCollectionTpl<PhaseSpec>::ImpulseDataVariant_t
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = ImpulseTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ImpulseDataBase<ImpulseDataTpl<PS, ImpulseCollectionTpl>, PS>;

        GALILEO_IMPULSE_DATA_TYPEDEF(Meta_t);

        using DataVariant_t = typename Collection_t::ImpulseDataVariant_t;

        DataVariant_t &toVariant() { return *static_cast<DataVariant_t *>(this); }
        const DataVariant_t &toVariant() const { return *static_cast<const DataVariant_t *>(this); }

        RobotData_t *robot() const { return galileo::impulse_robot_data(*this); }
        FrameIndex_t frame() const { return galileo::impulse_frame(*this); }
        ReferenceFrame_t type() const { return galileo::impulse_type_data(*this); }
        SE3_t jMf() const { return galileo::impulse_jMf(*this); }
        MatrixNcNv_t Jc() const { return galileo::impulse_Jc(*this); }
        Force_t f() const { return galileo::impulse_f(*this); }
        Force_t fext() const { return galileo::impulse_fext(*this); }
        MatrixNcNdx_t df_dx() const { return galileo::impulse_df_dx(*this); }
        MatrixNcNu_t df_du() const { return galileo::impulse_df_du(*this); }
        ActionMatrix_t fXj() const { return galileo::impulse_fXj(*this); }
        MatrixNcNv_t dv0_dq() const { return galileo::impulse_dv0_dq(*this); }
        MatrixNv_t dtau_dq() const { return galileo::impulse_dtau_dq(*this); }

        ImpulseDataTpl() : DataVariant_t() {}
        ImpulseDataTpl(const DataVariant_t &data_variant) : DataVariant_t(data_variant) {}
        template <typename DataDerived>
        ImpulseDataTpl(const ImpulseDataBase<DataDerived, PhaseSpec> &data)
            : Collection_t::ImpulseDataVariant_t((DataVariant_t) data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename DataVariant_t::types, DataDerived>) );
        }

        GENERIC_ACCESSOR(RobotData_t *, robot);
        GENERIC_ACCESSOR(FrameIndex_t, frame);
        GENERIC_ACCESSOR(ReferenceFrame_t, type);
        GENERIC_ACCESSOR(SE3_t, jMf);
        GENERIC_ACCESSOR(MatrixNcNv_t, Jc);
        GENERIC_ACCESSOR(Force_t, f);
        GENERIC_ACCESSOR(Force_t, fext);
        GENERIC_ACCESSOR(MatrixNcNdx_t, df_dx);
        GENERIC_ACCESSOR(MatrixNcNu_t, df_du);

        GENERIC_ACCESSOR(ActionMatrix_t, fXj);
        GENERIC_ACCESSOR(MatrixNcNv_t, dv0_dq);
        GENERIC_ACCESSOR(MatrixNv_t, dtau_dq);
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseModelTpl : public ImpulseModelBase<ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>, PhaseSpec>,
                             ImpulseCollectionTpl<PhaseSpec>::ImpulseModelVariant_t
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = ImpulseTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ImpulseModelBase<ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>, PS>;

        using ModelVariant_t = typename Collection_t::ImpulseModelVariant_t;

        using DimNC_t = typename traits<Meta_t>::DimNC_t;

        ImpulseModelTpl() : ModelVariant_t() {}
        ImpulseModelTpl(const ModelVariant_t &model_variant) : ModelVariant_t(model_variant) {}
        template <typename ModelDerived>
        ImpulseModelTpl(const ImpulseModelBase<ModelDerived, PhaseSpec> &model)
            : Base(model.get_id(), model.get_type(), DimNC_t(model.get_nc())),
              ModelVariant_t((ModelVariant_t) model.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename ModelVariant_t::types, ModelDerived>) );
        }

        ModelVariant_t &toVariant() { return *static_cast<ModelVariant_t *>(this); }
        const ModelVariant_t &toVariant() const { return *static_cast<const ModelVariant_t *>(this); }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            galileo::impulse_calc_zeroth_order(*this, data, x.derived());
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            galileo::impulse_calc_first_order(*this, data, x.derived());
        }

        Data_t createData(RobotData_t *const robot) const { return galileo::impulse_create_data(*this, robot); }

        template <typename ForceVectorType>
        void updateForce(Data_t &data, const Eigen::MatrixBase<ForceVectorType> &force) const
        {
            galileo::impulse_update_force(*this, data, force.derived());
        }

        using Base::setZeroForce;
        using Base::setZeroForceDiff;
        using Base::updateForceDiff;

        template <typename MatrixNcNdxType>
        void updateForceDiffImpl(Data_t &data, const Eigen::MatrixBase<MatrixNcNdxType> &df_dx) const
        {
            galileo::impulse_update_force_diff(*this, data, df_dx.derived());
        }

        void setZeroForceImpl(Data_t &data) const { galileo::impulse_set_zero_force(*this, data); }
        void setZeroForceDiffImpl(Data_t &data) const { galileo::impulse_set_zero_force_diff(*this, data); }
        FrameIndex_t get_id_impl() const { return galileo::impulse_get_id(*this); }
        void set_id_impl(const FrameIndex_t &id) { galileo::impulse_set_id(*this, id); }
        ReferenceFrame_t get_type_impl() const { return galileo::impulse_get_type(*this); }
        void set_type_impl(const ReferenceFrame_t &type) { galileo::impulse_set_type(*this, type); }
        int get_nc_impl() const { return galileo::impulse_get_nc(*this); }

        using Base::get_id;
        using Base::set_id;
        using Base::get_type;
        using Base::set_type;
        using Base::get_nc;

    }; // struct ImpulseModelTpl

} // namespace galileo

#endif // __galileo_multibody_impulses_impulse_generic_hpp__
