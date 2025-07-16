#ifndef __galileo_core_constraints_equality_constraint_generic_hpp__
#define __galileo_core_constraints_equality_constraint_generic_hpp__

#include "galileo/core/constraints/equality/constraint-base.hpp"
#include "galileo/core/constraints/equality/constraint-collection.hpp"
#include "galileo/core/constraints/equality/constraint-visitors.hxx"
#include "galileo/core/constraints/equality/fwd.hpp"

#include <boost/mpl/contains.hpp>

namespace galileo
{

    template <
        typename PhaseSpec,
        template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintTpl;

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = ConstraintCollectionTpl<PS>;
        using Model_t = ConstraintModelTpl<PS, ConstraintCollectionTpl>;
        using Data_t = ConstraintDataTpl<PS, ConstraintCollectionTpl>;

        using DimNH_t = DimensionTpl<Eigen::Dynamic>;
        static constexpr int NH = DimNH_t::Value;

        using H_t = Eigen::GMatrix<typename PS::VarScalar, NH, 1, PS::Options>;
        using Hx_t = Eigen::GMatrix<typename PS::VarScalar, NH, PS::DimNDX_t::Value, PS::Options>;
        using Hu_t = Eigen::GMatrix<typename PS::VarScalar, NH, PS::DimNU_t::Value, PS::Options>;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintDataTpl
        : public ConstraintDataBase<ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>, PhaseSpec>,
          ConstraintCollectionTpl<PhaseSpec>::ConstraintDataVariant_t
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ConstraintDataBase<ConstraintDataTpl<PS, ConstraintCollectionTpl>, PS>;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(Meta_t);

        using DataVariant_t = typename Collection_t::ConstraintDataVariant_t;

        DataVariant_t &toVariant()
        {
            return *static_cast<DataVariant_t *>(this);
        }
        const DataVariant_t &toVariant() const
        {
            return *static_cast<const DataVariant_t *>(this);
        }

        H_t H() const
        {
            return galileo::constraint_H(*this);
        }

        Hx_t Hx() const
        {
            return galileo::constraint_Hx(*this);
        }

        Hu_t Hu() const
        {
            return galileo::constraint_Hu(*this);
        }

        ConstraintDataTpl()
            : DataVariant_t()
        {
        }

        ConstraintDataTpl(const DataVariant_t &data_variant)
            : DataVariant_t(data_variant)
        {
        }

        template <typename DataDerived>
        ConstraintDataTpl(const ConstraintDataBase<DataDerived, PhaseSpec> &data)
            : Collection_t::ConstraintDataVariant_t((DataVariant_t)data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename DataVariant_t::types, DataDerived>));
        }

        GENERIC_ACCESSOR(H_t, H);
        GENERIC_ACCESSOR(Hx_t, Hx);
        GENERIC_ACCESSOR(Hu_t, Hu);

    }; // struct ConstraintDataTpl

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintModelTpl
    : public ConstraintModelBase<ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>, PhaseSpec>,
                                ConstraintCollectionTpl<PhaseSpec>::ConstraintModelVariant_t
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ConstraintModelBase<ConstraintModelTpl<PS, ConstraintCollectionTpl>, PS>;

        using ModelVariant_t = typename Collection_t::ConstraintModelVariant_t;

        ModelVariant_t &toVariant()
        {
            return *static_cast<ModelVariant_t *>(this);
        }

        const ModelVariant_t &toVariant() const
        {
            return *static_cast<const ModelVariant_t *>(this);
        }

        ConstraintModelTpl()
            : ModelVariant_t()
        {
        }

        ConstraintModelTpl(const ModelVariant_t &model_variant)
            : ModelVariant_t(model_variant)
        {
        }

        template <typename ModelDerived>
        ConstraintModelTpl(const ConstraintModelBase<ModelDerived, PhaseSpec> &model)
            : Collection_t::ConstraintModelVariant_t((ModelVariant_t)model.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename ModelVariant_t::types, ModelDerived>));
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return galileo::constraint_create_data(*this, collector);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            galileo::constraint_calc_zeroth_order(*this, data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            galileo::constraint_calc_zeroth_order(*this, data, x.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            galileo::constraint_calc_first_order(*this, data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            galileo::constraint_calc_first_order(*this, data, x.derived());
        }

        using Base::get_ps;

        using Base::get_nh;
        using Base::get_nh_dim;

    }; // struct ConstraintModelTpl

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_generic_hpp__
