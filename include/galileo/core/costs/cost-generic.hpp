#ifndef __galileo_core_costs_cost_generic_hpp__
#define __galileo_core_costs_cost_generic_hpp__

#include "galileo/core/costs/fwd.hpp"
#include "galileo/core/costs/cost-base.hpp"
#include "galileo/core/costs/cost-collection.hpp"
#include "galileo/core/costs/cost-basic-visitors.hxx"

#include <boost/mpl/contains.hpp>

namespace galileo
{

    namespace core
    {

        template <
            typename PhaseSpec,
            template <typename PS> class CostCollectionTpl>
        struct CostTpl;

        template <typename PhaseSpec,
                  template <typename PS> class CostCollectionTpl>
        struct traits<CostTpl<PhaseSpec, CostCollectionTpl>>
        {
            using PS = PhaseSpec;
            using CostCollection = CostCollectionTpl<PS>;

            using CostDataDerived = CostDataTpl<PS, CostCollectionTpl>;
            using CostModelDerived = CostModelTpl<PS, CostCollectionTpl>;

            using L_t = typename PS::VarScalar;
            using Lx_t = Eigen::Matrix<typename PS::VarScalar, PS::NDX, 1, PS::Options>;
            using Lu_t = Eigen::Matrix<typename PS::VarScalar, PS::NU, 1, PS::Options>;
            using Lxx_t = Eigen::Matrix<typename PS::VarScalar, PS::NDX, PS::NDX, PS::Options>;
            using Lxu_t = Eigen::Matrix<typename PS::VarScalar, PS::NDX, PS::NU, PS::Options>;
            using Luu_t = Eigen::Matrix<typename PS::VarScalar, PS::NU, PS::NU, PS::Options>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class CostCollectionTpl>
        struct traits<CostDataTpl<PhaseSpec, CostCollectionTpl>>
        {
            using PS = PhaseSpec;
            using CostDerived = CostTpl<PS, CostCollectionTpl>;
            using CostDataDerived = typename traits<CostDerived>::CostDataDerived;
            using CostModelDerived = typename traits<CostDerived>::CostModelDerived;
        };

        template <typename PhaseSpec,
                  template <typename PS> class CostCollectionTpl>
        struct traits<CostModelTpl<PhaseSpec, CostCollectionTpl>>
        {
            using PS = PhaseSpec;
            using CostDerived = CostTpl<PS, CostCollectionTpl>;
            using CostDataDerived = typename traits<CostDerived>::CostDataDerived;
            using CostModelDerived = typename traits<CostDerived>::CostModelDerived;
        };

        template <typename PhaseSpec,
                  template <typename PS> class CostCollectionTpl>
        struct CostDataTpl : public CostDataBase<CostDataTpl<PhaseSpec, CostCollectionTpl>, PhaseSpec>,
                             CostCollectionTpl<PhaseSpec>::CostDataVariant
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using CostDerived = CostTpl<PS, CostCollectionTpl>;
            using CostModelDerived = typename traits<CostDerived>::CostModelDerived;
            using CostDataDerived = typename traits<CostDerived>::CostDataDerived;

            GALILEO_COST_DATA_TYPEDEF(CostDerived);

            using CostCollection = CostCollectionTpl<PS>;
            using CostDataVariant = typename CostCollection::CostDataVariant;

            CostDataVariant &toVariant()
            {
                return *static_cast<CostDataVariant *>(this);
            }
            const CostDataVariant &toVariant() const
            {
                return *static_cast<const CostDataVariant *>(this);
            }

            L_t L() const
            {
                return galileo::core::cost_L(*this);
            }

            Lx_t Lx() const
            {
                return galileo::core::cost_Lx(*this);
            }

            Lu_t Lu() const
            {
                return galileo::core::cost_Lu(*this);
            }

            Lxx_t Lxx() const
            {
                return galileo::core::cost_Lxx(*this);
            }

            Lxu_t Lxu() const
            {
                return galileo::core::cost_Lxu(*this);
            }

            Luu_t Luu() const
            {
                return galileo::core::cost_Luu(*this);
            }

            CostDataTpl()
                : CostDataVariant()
            {
            }

            CostDataTpl(const CostDataVariant &cost_data_variant)
                : CostDataVariant(cost_data_variant)
            {
            }

            template <typename CostDataDerived>
            CostDataTpl(const CostDataBase<CostDataDerived, PhaseSpec> &cost_data)
                : CostCollection::CostDataVariant((CostDataVariant)cost_data.derived())
            {
                BOOST_MPL_ASSERT((boost::mpl::contains<typename CostDataVariant::types, CostDataDerived>));
            }

            GENERIC_ACCESSOR(L_t, L);
            GENERIC_ACCESSOR(Lx_t, Lx);
            GENERIC_ACCESSOR(Lu_t, Lu);
            GENERIC_ACCESSOR(Lxx_t, Lxx);
            GENERIC_ACCESSOR(Lxu_t, Lxu);
            GENERIC_ACCESSOR(Luu_t, Luu);

        }; // struct CostDataTpl

        template <typename PhaseSpec,
                  template <typename PS> class CostCollectionTpl>
        struct CostModelTpl : public CostModelBase<CostModelTpl<PhaseSpec, CostCollectionTpl>, PhaseSpec>,
                              CostCollectionTpl<PhaseSpec>::CostModelVariant
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using CostDerived = CostTpl<PS, CostCollectionTpl>;
            using CostModelDerived = typename traits<CostDerived>::CostModelDerived;
            using CostDataDerived = typename traits<CostDerived>::CostDataDerived;

            using CostCollection = CostCollectionTpl<PS>;
            using CostModelVariant = typename CostCollection::CostModelVariant;

            CostModelTpl()
                : CostModelVariant()
            {
            }

            CostModelTpl(const CostModelVariant &cost_model_variant)
                : CostModelVariant(cost_model_variant)
            {
            }

            template <typename CostModelDerived>
            CostModelTpl(const CostModelBase<CostModelDerived, PhaseSpec> &cost_model)
                : CostCollection::CostModelVariant((CostModelVariant)cost_model.derived())
            {
                BOOST_MPL_ASSERT((boost::mpl::contains<typename CostModelVariant::types, CostModelDerived>));
            }

            CostModelVariant &toVariant()
            {
                return *static_cast<CostModelVariant *>(this);
            }

            const CostModelVariant &toVariant() const
            {
                return *static_cast<const CostModelVariant *>(this);
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calc(CostDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                galileo::core::cost_calc_zeroth_order(*this, data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calc(CostDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
            {
                galileo::core::cost_calc_zeroth_order(*this, data, x.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(CostDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                galileo::core::cost_calc_first_order(*this, data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calcDiff(CostDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
            {
                galileo::core::cost_calc_first_order(*this, data, x.derived());
            }

            template <typename DataCollector>
            auto createData(DataCollector *const collector)
            {
                return galileo::core::cost_create_data(*this, collector);
            }

        }; // struct CostModelTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_generic_hpp__