#ifndef __galileo_predictive_phases_phase_generic_hpp__
#define __galileo_predictive_phases_phase_generic_hpp__

#include "galileo/predictive/phases/fwd.hpp"
#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-collection.hpp"
#include "galileo/predictive/phases/phase-visitors.hxx"

#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <
        typename PhaseSpec,
        template <typename PS> class PhaseCollectionTpl>
    struct PhaseTpl;

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    struct traits<PhaseTpl<PhaseSpec, PhaseCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseTpl<PS, PhaseCollectionTpl>;
        using Collection_t = PhaseCollectionTpl<PS>;
        using Model_t = PhaseModelTpl<PS, PhaseCollectionTpl>;
        using Data_t = PhaseDataTpl<PS, PhaseCollectionTpl>;
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    struct traits<PhaseDataTpl<PhaseSpec, PhaseCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseTpl<PS, PhaseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    struct traits<PhaseModelTpl<PhaseSpec, PhaseCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseTpl<PS, PhaseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    struct PhaseDataTpl : public PhaseDataBase<PhaseDataTpl<PhaseSpec, PhaseCollectionTpl>, PhaseSpec>,
                          PhaseCollectionTpl<PhaseSpec>::PhaseDataVariant_t
    {
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseTpl<PS, PhaseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseDataBase<PhaseDataTpl<PS, PhaseCollectionTpl>, PS>;

        using DataVariant_t = typename Collection_t::PhaseDataVariant_t;

        DataVariant_t &toVariant()
        {
            return *static_cast<DataVariant_t *>(this);
        }
        const DataVariant_t &toVariant() const
        {
            return *static_cast<const DataVariant_t *>(this);
        }

        PhaseDataTpl()
            : DataVariant_t()
        {
        }

        PhaseDataTpl(const DataVariant_t &data_variant)
            : DataVariant_t(data_variant)
        {
        }

        template <typename DataDerived>
        PhaseDataTpl(const PhaseDataBase<DataDerived, PhaseSpec> &data)
            : Collection_t::PhaseDataVariant_t((DataVariant_t)data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename DataVariant_t::types, DataDerived>));
        }

    }; // struct PhaseDataTpl

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    struct PhaseModelTpl : public PhaseModelBase<PhaseModelTpl<PhaseSpec, PhaseCollectionTpl>, PhaseSpec>,
                           PhaseCollectionTpl<PhaseSpec>::PhaseModelVariant_t
    {
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseTpl<PS, PhaseCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseModelBase<PhaseModelTpl<PS, PhaseCollectionTpl>, PS>;

        using ModelVariant_t = typename Collection_t::PhaseModelVariant_t;

        PhaseModelTpl()
            : ModelVariant_t()
        {
        }

        PhaseModelTpl(const ModelVariant_t &model_variant)
            : ModelVariant_t(model_variant)
        {
        }

        template <typename ModelDerived>
        PhaseModelTpl(const PhaseModelBase<ModelDerived, PhaseSpec> &model)
            : Base(model.get_ps()),
              ModelVariant_t((ModelVariant_t)model.derived())
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
                         const typename PS::NumScalar tol) const
        {
            galileo::phase_quasi_static(*this, data, xs.derived(), ws.derived(), maxiter, tol);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return galileo::phase_create_data(*this, collector);
        }

    }; // struct PhaseModelTpl

} // namespace galileo

#endif // __galileo_predictive_phases_phase_generic_hpp__
