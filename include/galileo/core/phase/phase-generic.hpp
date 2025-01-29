#ifndef __galileo_core_phase_generic_hpp__
#define __galileo_core_phase_generic_hpp__

#include "galileo/core/phase/phase-collection.hpp"
#include "galileo/core/phase/phase-basic-visitors.hxx"
#include "galileo/utils/aligned-vector.hpp"

#include <boost/mpl/contains.hpp>

namespace galileo
{

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class PhaseCollectionTpl = PhaseCollectionDefaultTpl>
    struct PhaseTpl;
    using Phase = PhaseTpl<context::Scalar>;

    template <typename _Scalar, int _Options, template <typename S, int O> class PhaseCollectionTpl>
    struct traits<PhaseTpl<_Scalar, _Options, PhaseCollectionTpl>>
    {
        enum
        {
            Options = _Options,
            NX = Eigen::Dynamic, // Dynamic because unknown at compile time
            NDX = Eigen::Dynamic,
            NU = Eigen::Dynamic,
            NH = Eigen::Dynamic,
            NG = Eigen::Dynamic
        };

        typedef _Scalar Scalar;
        typedef PhaseCollectionTpl<Scalar, Options> PhaseCollection;
        typedef PhaseDataTpl<Scalar, Options, PhaseCollectionTpl> PhaseDataDerived;
        typedef PhaseModelTpl<Scalar, Options, PhaseCollectionTpl> PhaseModelDerived;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class PhaseCollectionTpl>
    struct traits<PhaseDataTpl<_Scalar, _Options, PhaseCollectionTpl>>
    {
        using PhaseDerived = PhaseTpl<_Scalar, _Options, PhaseCollectionTpl>;
        using Scalar = typename traits<PhaseDerived>::Scalar;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class PhaseCollectionTpl>
    struct traits<PhaseModelTpl<_Scalar, _Options, PhaseCollectionTpl>>
    {
        using PhaseDerived = PhaseTpl<_Scalar, _Options, PhaseCollectionTpl>;
        using Scalar = typename traits<PhaseDerived>::Scalar;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class PhaseCollectionTpl>
    struct PhaseDataTpl
        : public PhaseDataBase<PhaseDataTpl<_Scalar, _Options, PhaseCollectionTpl>>,
          PhaseCollectionTpl<_Scalar, _Options>::PhaseDataVariant
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PhaseDerived = PhaseTpl<_Scalar, _Options, PhaseCollectionTpl>;
        using Base = PhaseDataBase<PhaseDataTpl>;

        GALILEO_PHASE_DATA_TYPEDEF_TEMPLATE(PhaseDerived);

        using PhaseCollection = PhaseCollectionTpl<_Scalar, _Options>;
        using PhaseDataVariant = typename PhaseCollection::PhaseDataVariant;

        using Base::operator==;
        using Base::operator!=;

        PhaseDataVariant &toVariant()
        {
            return *static_cast<PhaseDataVariant *>(this);
        }
        const PhaseDataVariant &toVariant() const
        {
            return *static_cast<const PhaseDataVariant *>(this);
        }

        PhaseDataTpl()
            : PhaseDataVariant()
        {
        }

        PhaseDataTpl(const PhaseDataVariant &jdata_variant)
            : PhaseDataVariant(jdata_variant)
        {
        }

        template <typename PhaseDataDerived>
        PhaseDataTpl(const PhaseDataBase<PhaseDataDerived> &jdata)
            : PhaseCollection::PhaseDataVariant((PhaseDataVariant)jdata.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename PhaseDataVariant::types, PhaseDataDerived>));
        }

        // Define all the standard accessors
    };

    template <
        typename NewScalar,
        typename Scalar,
        int Options,
        template <typename S, int O> class PhaseCollectionTpl>
    struct CastType<NewScalar, PhaseModelTpl<Scalar, Options, PhaseCollectionTpl>>
    {
        using type = PhaseModelTpl<NewScalar, Options, PhaseCollectionTpl>;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class PhaseCollectionTpl>
    struct PhaseModelTpl
        : PhaseModelBase<PhaseModelTpl<_Scalar, _Options, PhaseCollectionTpl>>,
          PhaseCollectionTpl<_Scalar, _Options>::PhaseModelVariant
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PhaseDerived = PhaseTpl<_Scalar, _Options, PhaseCollectionTpl>;

        GALILEO_PHASE_TYPEDEF_TEMPLATE(PhaseDerived);

        using PhaseCollection = PhaseCollectionTpl<Scalar, Options>;
        using PhaseDataVariant = typename PhaseCollection::PhaseDataVariant;
        using PhaseModelVariant = typename PhaseCollection::PhaseModelVariant;

        PhaseModelTpl()
            : PhaseModelVariant()
        {
        }

        PhaseModelTpl(const PhaseModelVariant &jmodel_variant)
            : PhaseCollection::PhaseModelVariant(jmodel_variant)
        {
        }

        template <typename PhaseModelDerived>
        PhaseModelTpl(const PhaseModelBase<PhaseModelDerived> &jmodel)
            : PhaseModelVariant((PhaseModelVariant)jmodel.derived())
        {
            BOOST_MPL_ASSERT(
                (boost::mpl::contains<typename PhaseModelVariant::types, PhaseModelDerived>));
        }

        PhaseModelVariant &toVariant()
        {
            return *static_cast<PhaseModelVariant *>(this);
        }

        const PhaseModelVariant &toVariant() const
        {
            return *static_cast<const PhaseModelVariant *>(this);
        }

        PhaseDataDerived createData() const
        {
            return ::galileo::createData<Scalar, Options, PhaseCollectionTpl>(*this);
        }

        /// \returns An expression of *this with the Scalar type casted to NewScalar.
        template <typename NewScalar>
        PhaseModelTpl<NewScalar, Options, PhaseCollectionTpl> cast() const
        {
            return cast_joint<NewScalar, Scalar, Options, PhaseCollectionTpl>(*this);
        }
    };

    using PhaseModelVector = typename GALILEO_ALIGNED_STD_VECTOR(PhaseData);
    using PhaseDataVector = typename GALILEO_ALIGNED_STD_VECTOR(PhaseModel);

} // namespace galileo

#endif // __galileo_core_phase_generic_hpp__
