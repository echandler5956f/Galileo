#ifndef __galileo_predictive_phases_phase_generic_hpp__
#define __galileo_predictive_phases_phase_generic_hpp__

#include "galileo/predictive/phases/fwd.hpp"
#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-collection.hpp"
#include "galileo/predictive/phases/phase-visitors.hxx"

#include "galileo/core/basic-spec.hpp"

namespace galileo
{

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseTpl;

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct traits<PhaseTpl<BasicSpec, PhaseCollectionTpl>>
    {
        using BS = BasicSpec;
        using SpecOfBaseClass = BS;

        GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BS);

        using Meta_t = PhaseTpl<BS, PhaseCollectionTpl>;
        using Collection_t = PhaseCollectionTpl<BS>;
        using Model_t = PhaseModelTpl<BS, PhaseCollectionTpl>;
        using Data_t = PhaseDataTpl<BS, PhaseCollectionTpl>;

        using XNext_t = VectorX_t;
        using XNextx_t = MatrixX_t;
        using XNextw_t = MatrixX_t;
        using L_t = VarScalar;
        using Lx_t = VectorX_t;
        using Lw_t = VectorX_t;
        using Lxx_t = MatrixX_t;
        using Lxw_t = MatrixX_t;
        using Lww_t = MatrixX_t;
        using H_t = VectorX_t;
        using Hx_t = MatrixX_t;
        using Hw_t = MatrixX_t;
        using G_t = VectorX_t;
        using Gx_t = MatrixX_t;
        using Gw_t = MatrixX_t;
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct traits<PhaseDataTpl<BasicSpec, PhaseCollectionTpl>>
    {
        using BS = BasicSpec;
        using SpecOfBaseClass = BS;

        using Meta_t = PhaseTpl<BS, PhaseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct traits<PhaseModelTpl<BasicSpec, PhaseCollectionTpl>>
    {
        using BS = BasicSpec;
        using SpecOfBaseClass = BS;

        using Meta_t = PhaseTpl<BS, PhaseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseDataTpl : public PhaseDataBase<PhaseDataTpl<BasicSpec, PhaseCollectionTpl>, BasicSpec>,
                          PhaseCollectionTpl<BasicSpec>::PhaseDataVariant_t
    {
        using BS = BasicSpec;

        GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BS);

        using Meta_t = PhaseTpl<BS, PhaseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseDataBase<PhaseDataTpl<BS, PhaseCollectionTpl>, BS>;

        GALILEO_PHASE_DATA_TYPEDEF(Meta_t);

        using DataVariant_t = typename Collection_t::PhaseDataVariant_t;

        DataVariant_t &toVariant() { return *static_cast<DataVariant_t *>(this); }
        const DataVariant_t &toVariant() const { return *static_cast<const DataVariant_t *>(this); }

        PhaseDataTpl() : DataVariant_t() {}
        PhaseDataTpl(const DataVariant_t &data_variant) : DataVariant_t(data_variant) {}
        template <typename DataDerived>
        PhaseDataTpl(const PhaseDataBase<DataDerived, BS> &data)
            : Collection_t::PhaseDataVariant_t((DataVariant_t) data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename DataVariant_t::types, DataDerived>) );
        }

        XNext_t &XNext_at_seg_i(const int i) const { return galileo::phase_XNext_at_seg_i(*this, i); }
        XNextx_t &XNextx_at_seg_i(const int i) const { return galileo::phase_XNextx_at_seg_i(*this, i); }
        XNextw_t &XNextw_at_seg_i(const int i) const { return galileo::phase_XNextw_at_seg_i(*this, i); }
        L_t &L_at_seg_i(const int i) const { return galileo::phase_L_at_seg_i(*this, i); }
        Lx_t &Lx_at_seg_i(const int i) const { return galileo::phase_Lx_at_seg_i(*this, i); }
        Lw_t &Lw_at_seg_i(const int i) const { return galileo::phase_Lw_at_seg_i(*this, i); }
        Lxx_t &Lxx_at_seg_i(const int i) const { return galileo::phase_Lxx_at_seg_i(*this, i); }
        Lxw_t &Lxw_at_seg_i(const int i) const { return galileo::phase_Lxw_at_seg_i(*this, i); }
        Lww_t &Lww_at_seg_i(const int i) const { return galileo::phase_Lww_at_seg_i(*this, i); }
        H_t &H_at_seg_i(const int i) const { return galileo::phase_H_at_seg_i(*this, i); }
        Hx_t &Hx_at_seg_i(const int i) const { return galileo::phase_Hx_at_seg_i(*this, i); }
        Hw_t &Hw_at_seg_i(const int i) const { return galileo::phase_Hw_at_seg_i(*this, i); }
        G_t &G_at_seg_i(const int i) const { return galileo::phase_G_at_seg_i(*this, i); }
        Gx_t &Gx_at_seg_i(const int i) const { return galileo::phase_Gx_at_seg_i(*this, i); }
        Gw_t &Gw_at_seg_i(const int i) const { return galileo::phase_Gw_at_seg_i(*this, i); }

        /////////////////////////////////////////////////////////////

        XNext_t &XNext_at_seg_i_accessor(const int i) { return XNext_at_seg_i(i); }
        const XNext_t &XNext_at_seg_i_accessor(const int i) const { return XNext_at_seg_i(i); }
        XNextx_t &XNextx_at_seg_i_accessor(const int i) { return XNextx_at_seg_i(i); }
        const XNextx_t &XNextx_at_seg_i_accessor(const int i) const { return XNextx_at_seg_i(i); }
        XNextw_t &XNextw_at_seg_i_accessor(const int i) { return XNextw_at_seg_i(i); }
        const XNextw_t &XNextw_at_seg_i_accessor(const int i) const { return XNextw_at_seg_i(i); }
        L_t &L_at_seg_i_accessor(const int i) { return L_at_seg_i(i); }
        const L_t &L_at_seg_i_accessor(const int i) const { return L_at_seg_i(i); }
        Lx_t &Lx_at_seg_i_accessor(const int i) { return Lx_at_seg_i(i); }
        const Lx_t &Lx_at_seg_i_accessor(const int i) const { return Lx_at_seg_i(i); }
        Lw_t &Lw_at_seg_i_accessor(const int i) { return Lw_at_seg_i(i); }
        const Lw_t &Lw_at_seg_i_accessor(const int i) const { return Lw_at_seg_i(i); }
        Lxx_t &Lxx_at_seg_i_accessor(const int i) { return Lxx_at_seg_i(i); }
        const Lxx_t &Lxx_at_seg_i_accessor(const int i) const { return Lxx_at_seg_i(i); }
        Lxw_t &Lxw_at_seg_i_accessor(const int i) { return Lxw_at_seg_i(i); }
        const Lxw_t &Lxw_at_seg_i_accessor(const int i) const { return Lxw_at_seg_i(i); }
        Lww_t &Lww_at_seg_i_accessor(const int i) { return Lww_at_seg_i(i); }
        const Lww_t &Lww_at_seg_i_accessor(const int i) const { return Lww_at_seg_i(i); }
        H_t &H_at_seg_i_accessor(const int i) { return H_at_seg_i(i); }
        const H_t &H_at_seg_i_accessor(const int i) const { return H_at_seg_i(i); }
        Hx_t &Hx_at_seg_i_accessor(const int i) { return Hx_at_seg_i(i); }
        const Hx_t &Hx_at_seg_i_accessor(const int i) const { return Hx_at_seg_i(i); }
        Hw_t &Hw_at_seg_i_accessor(const int i) { return Hw_at_seg_i(i); }
        const Hw_t &Hw_at_seg_i_accessor(const int i) const { return Hw_at_seg_i(i); }
        G_t &G_at_seg_i_accessor(const int i) { return G_at_seg_i(i); }
        const G_t &G_at_seg_i_accessor(const int i) const { return G_at_seg_i(i); }
        Gx_t &Gx_at_seg_i_accessor(const int i) { return Gx_at_seg_i(i); }
        const Gx_t &Gx_at_seg_i_accessor(const int i) const { return Gx_at_seg_i(i); }
        Gw_t &Gw_at_seg_i_accessor(const int i) { return Gw_at_seg_i(i); }
        const Gw_t &Gw_at_seg_i_accessor(const int i) const { return Gw_at_seg_i(i); }

    }; // struct PhaseDataTpl

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseModelTpl : public PhaseModelBase<PhaseModelTpl<BasicSpec, PhaseCollectionTpl>, BasicSpec>,
                           PhaseCollectionTpl<BasicSpec>::PhaseModelVariant_t
    {
        using BS = BasicSpec;

        GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BS);

        using Meta_t = PhaseTpl<BS, PhaseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseModelBase<PhaseModelTpl<BS, PhaseCollectionTpl>, BS>;

        using ModelVariant_t = typename Collection_t::PhaseModelVariant_t;

        PhaseModelTpl() : ModelVariant_t() {}
        PhaseModelTpl(const ModelVariant_t &model_variant) : ModelVariant_t(model_variant) {}
        template <typename ModelDerived>
        PhaseModelTpl(const PhaseModelBase<ModelDerived, BS> &model)
            : Base(), ModelVariant_t((ModelVariant_t) model.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename ModelVariant_t::types, ModelDerived>) );
        }

        ModelVariant_t &toVariant() { return *static_cast<ModelVariant_t *>(this); }
        const ModelVariant_t &toVariant() const { return *static_cast<const ModelVariant_t *>(this); }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs,
                  const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            galileo::phase_calc_zeroth_order(*this, data, xs.derived(), ws.derived());
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            galileo::phase_calc_first_order(*this, data, xs.derived(), ws.derived());
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void quasiStatic(Data_t &data,
                         const Eigen::MatrixBase<StateMatrixType> &xs,
                         Eigen::MatrixBase<ControlParamMatrixType> &ws,
                         const int maxiter,
                         const NumScalar tol) const
        {
            galileo::phase_quasi_static(*this, data, xs.derived(), ws.derived(), maxiter, tol);
        }

        Data_t createData() const { return galileo::phase_create_data(*this); }

    }; // struct PhaseModelTpl

} // namespace galileo

#endif // __galileo_predictive_phases_phase_generic_hpp__
