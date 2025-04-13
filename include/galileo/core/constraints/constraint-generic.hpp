#ifndef __galileo_core_constraints_constraint_generic_hpp__
#define __galileo_core_constraints_constraint_generic_hpp__

#include "galileo/core/constraints/fwd.hpp"
#include "galileo/core/constraints/constraint-base.hpp"
#include "galileo/core/constraints/constraint-collection.hpp"
#include "galileo/core/constraints/constraint-basic-visitors.hxx"

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

        static constexpr int NH = Eigen::Dynamic;
        static constexpr int NG = Eigen::Dynamic;

        using H_t = Eigen::Matrix<typename PS::VarScalar, NH, 1, PS::Options>;
        using Hx_t = Eigen::Matrix<typename PS::VarScalar, NH, PS::NDX, PS::Options>;
        using Hu_t = Eigen::Matrix<typename PS::VarScalar, NH, PS::NU, PS::Options>;
        using G_t = Eigen::Matrix<typename PS::VarScalar, NG, 1, PS::Options>;
        using Gx_t = Eigen::Matrix<typename PS::VarScalar, NG, PS::NDX, PS::Options>;
        using Gu_t = Eigen::Matrix<typename PS::VarScalar, NG, PS::NU, PS::Options>;

        using BoundVector_t = Eigen::Matrix<typename PS::NumScalar, NG, 1, PS::Options>;
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
    struct ConstraintDataTpl : public ConstraintDataBase<ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>, PhaseSpec>,
                               ConstraintCollectionTpl<PhaseSpec>::DataVariant_t
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(Meta_t);

        using DataVariant_t = typename Collection_t::DataVariant_t;

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

        G_t G() const
        {
            return galileo::constraint_G(*this);
        }

        Gx_t Gx() const
        {
            return galileo::constraint_Gx(*this);
        }

        Gu_t Gu() const
        {
            return galileo::constraint_Gu(*this);
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
            : Collection_t::DataVariant_t((DataVariant_t)data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename DataVariant_t::types, DataDerived>));
        }

        GENERIC_ACCESSOR(H_t, H);
        GENERIC_ACCESSOR(Hx_t, Hx);
        GENERIC_ACCESSOR(Hu_t, Hu);
        GENERIC_ACCESSOR(G_t, G);
        GENERIC_ACCESSOR(Gx_t, Gx);
        GENERIC_ACCESSOR(Gu_t, Gu);

    }; // struct ConstraintDataTpl

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintModelTpl : public ConstraintModelBase<ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>, PhaseSpec>,
                                ConstraintCollectionTpl<PhaseSpec>::ModelVariant_t
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using ModelVariant_t = typename Collection_t::ModelVariant_t;

        using BoundVector_t = typename traits<Meta_t>::BoundVector_t;

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
            : Collection_t::ModelVariant_t((ModelVariant_t)model.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename ModelVariant_t::types, ModelDerived>));
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector)
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

        template <typename LowerBoundType, typename UpperBoundType>
        void updateBounds(const Eigen::MatrixBase<LowerBoundType> &lb,
                          const Eigen::MatrixBase<UpperBoundType> &ub)
        {
            galileo::constraint_update_bounds(*this, lb.derived(), ub.derived());
        }

        const BoundVector_t &lb() const
        {
            return galileo::constraint_lb(*this);
        }

        const BoundVector_t &ub() const
        {
            return galileo::constraint_ub(*this);
        }

        int ng_impl() const
        {
            return galileo::constraint_ng(*this);
        }

        int nh_impl() const
        {
            return galileo::constraint_nh(*this);
        }

    }; // struct ConstraintModelTpl

} // namespace galileo

#endif // __galileo_core_constraints_constraint_generic_hpp__