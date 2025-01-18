#ifndef __galileo_core_segment_generic_hpp__
#define __galileo_core_segment_generic_hpp__

#include "galileo/core/segment/segment-collection.hpp"
#include "galileo/core/segment/segment-basic-visitors.hxx"
#include "galileo/utils/aligned-vector.hpp"

#include <boost/mpl/contains.hpp>

namespace galileo
{

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class SegmentCollectionTpl = SegmentCollectionDefaultTpl>
    struct SegmentTpl;
    using Segment = SegmentTpl<context::Scalar>;

    template <typename _Scalar, int _Options, template <typename S, int O> class SegmentCollectionTpl>
    struct traits<SegmentTpl<_Scalar, _Options, SegmentCollectionTpl>>
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
        typedef SegmentCollectionTpl<Scalar, Options> SegmentCollection;
        typedef SegmentDataTpl<Scalar, Options, SegmentCollectionTpl> SegmentDataDerived;
        typedef SegmentModelTpl<Scalar, Options, SegmentCollectionTpl> SegmentModelDerived;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class SegmentCollectionTpl>
    struct traits<SegmentDataTpl<_Scalar, _Options, SegmentCollectionTpl>>
    {
        using SegmentDerived = SegmentTpl<_Scalar, _Options, SegmentCollectionTpl>;
        using Scalar = typename traits<SegmentDerived>::Scalar;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class SegmentCollectionTpl>
    struct traits<SegmentModelTpl<_Scalar, _Options, SegmentCollectionTpl>>
    {
        using SegmentDerived = SegmentTpl<_Scalar, _Options, SegmentCollectionTpl>;
        using Scalar = typename traits<SegmentDerived>::Scalar;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class SegmentCollectionTpl>
    struct SegmentDataTpl
        : public SegmentDataBase<SegmentDataTpl<_Scalar, _Options, SegmentCollectionTpl>>,
          SegmentCollectionTpl<_Scalar, _Options>::SegmentDataVariant
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using SegmentDerived = SegmentTpl<_Scalar, _Options, SegmentCollectionTpl>;
        using Base = SegmentDataBase<SegmentDataTpl>;

        GALILEO_SEGMENT_DATA_TYPEDEF_TEMPLATE(SegmentDerived);

        using SegmentCollection = SegmentCollectionTpl<_Scalar, _Options>;
        using SegmentDataVariant = typename SegmentCollection::SegmentDataVariant;

        using Base::operator==;
        using Base::operator!=;

        SegmentDataVariant &toVariant()
        {
            return *static_cast<SegmentDataVariant *>(this);
        }
        const SegmentDataVariant &toVariant() const
        {
            return *static_cast<const SegmentDataVariant *>(this);
        }

        SegmentDataTpl()
            : SegmentDataVariant()
        {
        }

        SegmentDataTpl(const SegmentDataVariant &jdata_variant)
            : SegmentDataVariant(jdata_variant)
        {
        }

        template <typename SegmentDataDerived>
        SegmentDataTpl(const SegmentDataBase<SegmentDataDerived> &jdata)
            : SegmentCollection::SegmentDataVariant((SegmentDataVariant)jdata.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename SegmentDataVariant::types, SegmentDataDerived>));
        }

        // Define all the standard accessors
    };

    template <
        typename NewScalar,
        typename Scalar,
        int Options,
        template <typename S, int O> class SegmentCollectionTpl>
    struct CastType<NewScalar, SegmentModelTpl<Scalar, Options, SegmentCollectionTpl>>
    {
        using type = SegmentModelTpl<NewScalar, Options, SegmentCollectionTpl>;
    };

    template <typename _Scalar, int _Options, template <typename S, int O> class SegmentCollectionTpl>
    struct SegmentModelTpl
        : SegmentModelBase<SegmentModelTpl<_Scalar, _Options, SegmentCollectionTpl>>,
          SegmentCollectionTpl<_Scalar, _Options>::SegmentModelVariant
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using SegmentDerived = SegmentTpl<_Scalar, _Options, SegmentCollectionTpl>;

        GALILEO_SEGMENT_TYPEDEF_TEMPLATE(SegmentDerived);

        using SegmentCollection = SegmentCollectionTpl<Scalar, Options>;
        using SegmentDataVariant = typename SegmentCollection::SegmentDataVariant;
        using SegmentModelVariant = typename SegmentCollection::SegmentModelVariant;

        SegmentModelTpl()
            : SegmentModelVariant()
        {
        }

        SegmentModelTpl(const SegmentModelVariant &jmodel_variant)
            : SegmentCollection::SegmentModelVariant(jmodel_variant)
        {
        }

        template <typename SegmentModelDerived>
        SegmentModelTpl(const SegmentModelBase<SegmentModelDerived> &jmodel)
            : SegmentModelVariant((SegmentModelVariant)jmodel.derived())
        {
            BOOST_MPL_ASSERT(
                (boost::mpl::contains<typename SegmentModelVariant::types, SegmentModelDerived>));
        }

        SegmentModelVariant &toVariant()
        {
            return *static_cast<SegmentModelVariant *>(this);
        }

        const SegmentModelVariant &toVariant() const
        {
            return *static_cast<const SegmentModelVariant *>(this);
        }

        SegmentDataDerived createData() const
        {
            return ::galileo::createData<Scalar, Options, SegmentCollectionTpl>(*this);
        }

        /// \returns An expression of *this with the Scalar type casted to NewScalar.
        template <typename NewScalar>
        SegmentModelTpl<NewScalar, Options, SegmentCollectionTpl> cast() const
        {
            return cast_joint<NewScalar, Scalar, Options, SegmentCollectionTpl>(*this);
        }
    };

    using SegmentModelVector = typename GALILEO_ALIGNED_STD_VECTOR(SegmentData);
    using SegmentDataVector = typename GALILEO_ALIGNED_STD_VECTOR(SegmentModel);

} // namespace galileo

#endif // __galileo_core_segment_generic_hpp__
