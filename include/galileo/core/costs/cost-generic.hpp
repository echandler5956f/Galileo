#ifndef __galileo_core_costs_cost_generic_hpp__
#define __galileo_core_costs_cost_generic_hpp__

#include "galileo/core/costs/cost-base.hpp"
#include "galileo/core/costs/cost-collection.hpp"
#include "galileo/core/costs/cost-visitors.hxx"

#include <boost/mpl/contains.hpp>

namespace galileo
{

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostTpl;

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct traits<CostTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = CostTpl<PS, CostCollectionTpl>;
        using Collection_t = CostCollectionTpl<PS>;
        using Model_t = CostModelTpl<PS, CostCollectionTpl>;
        using Data_t = CostDataTpl<PS, CostCollectionTpl>;

        using L_t = typename PS::VarScalar;
        using Lx_t = Eigen::GMatrix<typename PS::VarScalar, PS::NDX, 1, PS::Options>;
        using Lu_t = Eigen::GMatrix<typename PS::VarScalar, PS::NU, 1, PS::Options>;
        using Lxx_t = Eigen::GMatrix<typename PS::VarScalar, PS::NDX, PS::NDX, PS::Options>;
        using Lxu_t = Eigen::GMatrix<typename PS::VarScalar, PS::NDX, PS::NU, PS::Options>;
        using Luu_t = Eigen::GMatrix<typename PS::VarScalar, PS::NU, PS::NU, PS::Options>;
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct traits<CostDataTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = CostTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct traits<CostModelTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = CostTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostDataTpl : public CostDataBase<CostDataTpl<PhaseSpec, CostCollectionTpl>, PhaseSpec>,
                         CostCollectionTpl<PhaseSpec>::CostDataVariant_t
    {
        using PS = PhaseSpec;

        using Meta_t = CostTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = CostDataBase<CostDataTpl<PS, CostCollectionTpl>, PS>;

        GALILEO_COST_DATA_TYPEDEF(Meta_t);

        using DataVariant_t = typename Collection_t::CostDataVariant_t;

        DataVariant_t &toVariant() { return *static_cast<DataVariant_t *>(this); }
        const DataVariant_t &toVariant() const { return *static_cast<const DataVariant_t *>(this); }

        L_t L() const { return galileo::cost_L(*this); }
        Lx_t Lx() const { return galileo::cost_Lx(*this); }
        Lu_t Lu() const { return galileo::cost_Lu(*this); }
        Lxx_t Lxx() const { return galileo::cost_Lxx(*this); }
        Lxu_t Lxu() const { return galileo::cost_Lxu(*this); }
        Luu_t Luu() const { return galileo::cost_Luu(*this); }

        CostDataTpl() : DataVariant_t() {}
        CostDataTpl(const DataVariant_t &data_variant) : DataVariant_t(data_variant) {}
        template <typename CostDataType>
        CostDataTpl(const CostDataBase<CostDataType, PhaseSpec> &data)
            : Collection_t::CostDataVariant_t((DataVariant_t) data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename DataVariant_t::types, CostDataType>) );
        }

        GENERIC_ACCESSOR(L_t, L);
        GENERIC_ACCESSOR(Lx_t, Lx);
        GENERIC_ACCESSOR(Lu_t, Lu);
        GENERIC_ACCESSOR(Lxx_t, Lxx);
        GENERIC_ACCESSOR(Lxu_t, Lxu);
        GENERIC_ACCESSOR(Luu_t, Luu);

    }; // struct CostDataTpl

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostModelTpl : public CostModelBase<CostModelTpl<PhaseSpec, CostCollectionTpl>, PhaseSpec>,
                          CostCollectionTpl<PhaseSpec>::CostModelVariant_t
    {
        using PS = PhaseSpec;

        using Meta_t = CostTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = CostModelBase<CostModelTpl<PS, CostCollectionTpl>, PS>;

        using ModelVariant_t = typename Collection_t::CostModelVariant_t;

        CostModelTpl() : ModelVariant_t() {}
        CostModelTpl(const ModelVariant_t &model_variant) : ModelVariant_t(model_variant) {}
        template <typename CostModelType>
        CostModelTpl(const CostModelBase<CostModelType, PhaseSpec> &model)
            : Base(), Collection_t::CostModelVariant_t((ModelVariant_t) model.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename ModelVariant_t::types, CostModelType>) );
        }

        ModelVariant_t &toVariant() { return *static_cast<ModelVariant_t *>(this); }
        const ModelVariant_t &toVariant() const { return *static_cast<const ModelVariant_t *>(this); }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            galileo::cost_calc_zeroth_order(*this, data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            galileo::cost_calc_zeroth_order(*this, data, x.derived(), Blank());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            galileo::cost_calc_first_order(*this, data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            galileo::cost_calc_first_order(*this, data, x.derived(), Blank());
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return galileo::cost_create_data(*this, collector);
        }

    }; // struct CostModelTpl

} // namespace galileo

#endif // __galileo_core_costs_cost_generic_hpp__
