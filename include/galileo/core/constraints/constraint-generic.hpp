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
            typename VarScalar,
            typename NumScalar,
            int Options,
            template <typename V, typename N, int O> class ConstraintCollectionTpl>
        struct ConstraintTpl
        {
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class ConstraintCollectionTpl>
        struct traits<ConstraintTpl<_VarScalar, _NumScalar, _Options, ConstraintCollectionTpl>>
        {
            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NX = ConstraintCollectionTpl<VarScalar, NumScalar, Options>::NX;
            static constexpr int NU = ConstraintCollectionTpl<VarScalar, NumScalar, Options>::NU;
            static constexpr int NH = Eigen::Dynamic;
            static constexpr int NG = Eigen::Dynamic;

            using ConstraintDataDerived = ConstraintDataTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl>;
            using ConstraintModelDerived = ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl>;

            using H_t = Eigen::Matrix<VarScalar, NH, NX, Options>;
            using Hx_t = Eigen::Matrix<VarScalar, NH, NX, Options>;
            using Hu_t = Eigen::Matrix<VarScalar, NH, NU, Options>;
            using G_t = Eigen::Matrix<VarScalar, NG, NX, Options>;
            using Gx_t = Eigen::Matrix<VarScalar, NG, NX, Options>;
            using Gu_t = Eigen::Matrix<VarScalar, NG, NU, Options>;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class ConstraintCollectionTpl>
        struct traits<ConstraintDataTpl<_VarScalar, _NumScalar, _Options, ConstraintCollectionTpl>>
        {
            typedef ConstraintTpl<_VarScalar, _NumScalar, _Options, ConstraintCollectionTpl> ConstraintDerived;
            typedef typename traits<ConstraintDerived>::VarScalar VarScalar;
            typedef typename traits<ConstraintDerived>::NumScalar NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class ConstraintCollectionTpl>
        struct traits<ConstraintModelTpl<_VarScalar, _NumScalar, _Options, ConstraintCollectionTpl>>
        {
            typedef ConstraintTpl<_VarScalar, _NumScalar, _Options, ConstraintCollectionTpl> ConstraintDerived;
            typedef typename traits<ConstraintDerived>::VarScalar VarScalar;
            typedef typename traits<ConstraintDerived>::NumScalar NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class ConstraintCollectionTpl>
        struct ConstraintDataTpl : ConstraintDataBase<ConstraintDataTpl<_VarScalar, _NumScalar, _Options, ConstraintCollectionTpl>>,
                                   ConstraintCollectionTpl<VarScalar, NumScalar, Options>::ConstraintDataVariant
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ConstraintDerived = ConstraintTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl>;
            using Base = ConstraintDataBase<ConstraintDataTpl>;

            GALILEO_CONSTRAINT_DATA_TYPEDEF(ConstraintDerived);

            using ConstraintCollection = ConstraintCollectionTpl<VarScalar, NumScalar, Options>;
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

            H_t H_accessor()
            {
                return H();
            }

            Hx_t Hx_accessor()
            {
                return Hx();
            }

            Hu_t Hu_accessor()
            {
                return Hu();
            }

            G_t G_accessor()
            {
                return G();
            }

            Gx_t Gx_accessor()
            {
                return Gx();
            }

            Gu_t Gu_accessor()
            {
                return Gu();
            }

        }; // struct ConstraintDataTpl

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class ConstraintCollectionTpl>
        struct ConstraintModelTpl : ConstraintModelBase<ConstraintModelTpl<_VarScalar, _NumScalar, _Options, ConstraintCollectionTpl>>,
                                    ConstraintCollectionTpl<VarScalar, NumScalar, Options>::ConstraintModelVariant
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ConstraintDerived = ConstraintTpl<VarScalar, NumScalar, Options, ConstraintCollectionTpl>;
            using Base = ConstraintModelBase<ConstraintModelTpl>;

            GALILEO_CONSTRAINT_MODEL_TYPEDEF(ConstraintDerived);

            using ConstraintCollection = ConstraintCollectionTpl<VarScalar, NumScalar, Options>;
            using ConstraintDataVariant = typename ConstraintCollection::ConstraintDataVariant;
            using ConstraintModelVariant = typename ConstraintCollection::ConstraintModelVariant;

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

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                galileo::core::constraint_calc_first_order(*this, data, x.derived(), u.derived());
            }

        }; // struct ConstraintModelTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_generic_hpp__