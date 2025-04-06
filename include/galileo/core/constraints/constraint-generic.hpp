#ifndef __galileo_core_constraints_constraint_generic_hpp__
#define __galileo_core_constraints_constraint_generic_hpp__

#include "galileo/core/constraints/fwd.hpp"
#include "galileo/core/constraints/constraint-base.hpp"
#include "galileo/core/constraints/constraint-collection.hpp"

#include <boost/mpl/contains.hpp>

namespace galileo
{

    namespace core
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
            using ConstraintCollection = ConstraintCollectionTpl<PS>;

            static constexpr int NH = Eigen::Dynamic;
            static constexpr int NG = Eigen::Dynamic;

            using ConstraintDataDerived = ConstraintDataTpl<PS, ConstraintCollectionTpl>;
            using ConstraintModelDerived = ConstraintModelTpl<PS, ConstraintCollectionTpl>;

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
            using ConstraintDerived = ConstraintTpl<PS, ConstraintCollectionTpl>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ConstraintCollectionTpl>
        struct traits<ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>>
        {
            using PS = PhaseSpec;
            using ConstraintDerived = ConstraintTpl<PS, ConstraintCollectionTpl>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ConstraintCollectionTpl>
        struct ConstraintDataTpl : public ConstraintDataBase<ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>, PhaseSpec>,
                                   ConstraintCollectionTpl<PhaseSpec>::ConstraintDataVariant
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ConstraintDerived = ConstraintTpl<PS, ConstraintCollectionTpl>;
            using ConstraintModelDerived = typename traits<ConstraintDerived>::ConstraintModelDerived;
            using ConstraintDataDerived = typename traits<ConstraintDerived>::ConstraintDataDerived;

            GALILEO_CONSTRAINT_DATA_TYPEDEF(ConstraintDerived);

            using ConstraintCollection = ConstraintCollectionTpl<PS>;
            using ConstraintDataVariant = typename ConstraintCollection::ConstraintDataVariant;

            ConstraintDataVariant &toVariant()
            {
                return *static_cast<ConstraintDataVariant *>(this);
            }
            const ConstraintDataVariant &toVariant() const
            {
                return *static_cast<const ConstraintDataVariant *>(this);
            }

            H_t H() const
            {
                return galileo::core::constraint_H(*this);
            }

            Hx_t Hx() const
            {
                return galileo::core::constraint_Hx(*this);
            }

            Hu_t Hu() const
            {
                return galileo::core::constraint_Hu(*this);
            }

            G_t G() const
            {
                return galileo::core::constraint_G(*this);
            }

            Gx_t Gx() const
            {
                return galileo::core::constraint_Gx(*this);
            }

            Gu_t Gu() const
            {
                return galileo::core::constraint_Gu(*this);
            }

            ConstraintDataTpl()
                : ConstraintDataVariant()
            {
            }

            ConstraintDataTpl(const ConstraintDataVariant &constraint_data_variant)
                : ConstraintDataVariant(constraint_data_variant)
            {
            }

            template <typename ConstraintDataDerived>
            ConstraintDataTpl(const ConstraintDataBase<ConstraintDataDerived> &constraint_data)
                : ConstraintCollection::ConstraintDataVariant((ConstraintDataVariant)constraint_data.derived())
            {
                BOOST_MPL_ASSERT((boost::mpl::contains<typename ConstraintDataVariant::types, ConstraintDataDerived>));
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
        struct ConstraintModelTpl : ConstraintModelBase<ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>>,
                                    ConstraintCollectionTpl<PS>::ConstraintModelVariant
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ConstraintDerived = ConstraintTpl<PS, ConstraintCollectionTpl>;
            using ConstraintModelDerived = typename traits<ConstraintDerived>::ConstraintModelDerived;
            using ConstraintDataDerived = typename traits<ConstraintDerived>::ConstraintDataDerived;

            using ConstraintCollection = ConstraintCollectionTpl<PS>;
            using ConstraintModelVariant = typename ConstraintCollection::ConstraintModelVariant;

            using BoundVector_t = typename traits<ConstraintDerived>::BoundVector_t;

            ConstraintModelTpl()
                : ConstraintModelVariant()
            {
            }

            ConstraintModelTpl(const ConstraintModelVariant &constraint_model_variant)
                : ConstraintModelVariant(constraint_model_variant)
            {
            }

            template <typename ConstraintModelDerived>
            ConstraintModelTpl(const ConstraintModelBase<ConstraintModelDerived> &constraint_model)
                : ConstraintCollection::ConstraintModelVariant((ConstraintModelVariant)constraint_model.derived())
            {
                BOOST_MPL_ASSERT((boost::mpl::contains<typename ConstraintModelVariant::types, ConstraintModelDerived>));
            }

            ConstraintModelVariant &toVariant()
            {
                return *static_cast<ConstraintModelVariant *>(this);
            }

            const ConstraintModelVariant &toVariant() const
            {
                return *static_cast<const ConstraintModelVariant *>(this);
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calc(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                galileo::core::constraint_calc_zeroth_order(*this, data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calc(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
            {
                galileo::core::constraint_calc_zeroth_order(*this, data, x.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                galileo::core::constraint_calc_first_order(*this, data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
            {
                galileo::core::constraint_calc_first_order(*this, data, x.derived());
            }

            template <typename DataCollector>
            auto createData(DataCollector *const collector)
            {
                return galileo::core::constraint_create_data(*this, collector);
            }

            template <typename LowerBoundType, typename UpperBoundType>
            void updateBounds(const Eigen::MatrixBase<LowerBoundType> &lb,
                              const Eigen::MatrixBase<UpperBoundType> &ub)
            {
                galileo::core::constraint_update_bounds(*this, lb.derived(), ub.derived());
            }

            const BoundVector_t &lb() const
            {
                return galileo::core::constraint_lb(*this);
            }

            const BoundVector_t &ub() const
            {
                return galileo::core::constraint_ub(*this);
            }

            int ng_impl() const
            {
                return galileo::core::constraint_ng(*this);
            }

            int nh_impl() const
            {
                return galileo::core::constraint_nh(*this);
            }

        }; // struct ConstraintModelTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_generic_hpp__