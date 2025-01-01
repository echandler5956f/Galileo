#ifndef __galileo_core_segment_model_base_hpp__
#define __galileo_core_segment_model_base_hpp__

#include "galileo/core/segment/segment-base.hpp"
#include "galileo/utils/aligned-vector.hpp"
#include <limits>

#define GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, TYPENAME)               \
    typedef TYPENAME traits<Segment>::Scalar Scalar;                           \
    typedef TYPENAME traits<Segment>::SegmentModelDerived SegmentModelDerived; \
    typedef TYPENAME traits<Segment>::SegmentDataDerived SegmentDataDerived;   \
    typedef TYPENAME traits<Segment>::NodeModel_t NodeModel_t;                 \
    typedef TYPENAME traits<Segment>::NodeData_t NodeData_t;                   \
    typedef TYPENAME traits<Segment>::State_t State_t;                         \
    typedef TYPENAME traits<Segment>::Control_t Control_t;                     \
    enum                                                                       \
    {                                                                          \
        Options = traits<Segment>::Options,                                    \
        NX = traits<Segment>::NX,                                              \
        NDX = traits<Segment>::NDX,                                            \
        NU = traits<Segment>::NU,                                              \
        NH = traits<Segment>::NH,                                              \
        NG = traits<Segment>::NG,                                              \
        NC = traits<Segment>::NC                                               \
    };                                                                         \
    typedef TYPENAME traits<Segment>::X_t X_t;                                 \
    typedef TYPENAME traits<Segment>::dX_t dX_t;                               \
    typedef TYPENAME traits<Segment>::U_t U_t;                                 \
    typedef TYPENAME traits<Segment>::H_t H_t;                                 \
    typedef TYPENAME traits<Segment>::G_t G_t;                                 \
    typedef TYPENAME traits<Segment>::C_t C_t;

#ifdef __clang__

#define GALILEO_SEGMENT_TYPEDEF(Segment) \
    GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, GALILEO_EMPTY_ARG)
#define GALILEO_SEGMENT_TYPEDEF_TEMPLATE(Segment) \
    GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, typename)

#elif (__GNUC__ == 4) && (__GNUC_MINOR__ == 4) && (__GNUC_PATCHLEVEL__ == 2)

#define GALILEO_SEGMENT_TYPEDEF(Segment) \
    GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, GALILEO_EMPTY_ARG)
#define GALILEO_SEGMENT_TYPEDEF_TEMPLATE(Segment) \
    GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, typename)

#else

#define GALILEO_SEGMENT_TYPEDEF(Segment) GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, typename)
#define GALILEO_SEGMENT_TYPEDEF_TEMPLATE(Segment) \
    GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, typename)

#endif

#define GALILEO_SEGMENT_CAST_TYPE_SPECIALIZATION(SegmentModelTpl) \
    template <typename Scalar, typename NewScalar>                \
    struct CastType<NewScalar, SegmentModelTpl<Scalar>>           \
    {                                                             \
        typedef SegmentModelTpl<NewScalar> type;                  \
    }

namespace galileo
{
    template <typename Derived>
    class SegmentModelBase : NumericalBase<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::SegmentDerived SegmentDerived;
        GALILEO_SEGMENT_TYPEDEF_TEMPLATE(SegmentDerived);

        SegmentModelDerived &derived()
        {
            return *static_cast<Derived *>(this);
        }

        const SegmentModelDerived &derived() const
        {
            return *static_cast<const Derived *>(this);
        }

        SegmentDataDerived createData() const
        {
            return derived().createData();
        }

        template <typename StateMatrixType, typename ControlMatrixType>
        void calc(SegmentDataDerived &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs,
                  const Eigen::MatrixBase<ControlMatrixType> &us) const
        {
            derived().calc(data, xs.derived(), us.derived());
        }

        template <typename StateMatrixType, typename ControlMatrixType>
        void calc(SegmentDataDerived &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs) const
        {
            derived().calc(data, xs.derived());
        }

        template <typename StateMatrixType, typename ControlMatrixType>
        void calcDiff(SegmentDataDerived &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlMatrixType> &us) const
        {
            derived().calcDiff(data, xs.derived(), us.derived());
        }

        template <typename StateMatrixType, typename ControlMatrixType>
        void calcDiff(SegmentDataDerived &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs) const
        {
            derived().calcDiff(data, xs.derived());
        }

        template <typename NewScalar>
        typename CastType<NewScalar, Derived>::type cast() const
        {
            return derived().template cast<NewScalar>();
        }

    protected:
        // Default constructor: protected.
        // Prevent the construction of stand-alone SegmentModelBase.
        inline SegmentModelBase() : period_(std::numeric_limits<double>::quiet_NaN), num_nodes_(std::numeric_limits<size_t>::max()), integ_constraint_size_(std::numeric_limits<size_t>::max())
        {
        }

        // Copy constructor: protected.
        // Copy of stand-alone SegmentModelBase are prevented, but can be used from inheriting
        // objects. Copy is done by calling copy operator.
        inline SegmentModelBase(const SegmentModelBase &clone)
        {
            *this = clone;
        }

        // Copy operator: protected.
        // Copy of stand-alone SegmentModelBase are prevented, but can be used from inheriting
        // objects.
        inline SegmentModelBase &operator=(const SegmentModelBase &clone)
        {
            node_models_ = clone.node_models_;
            state_desc_ = clone.state_desc_;
            control_desc_ = clone.control_desc_;
            period_ = clone.period_;
            num_nodes_ = clone.num_nodes_;
            nc_ = clone.nc_;
            return *this;
        }

        GALILEO_ALIGNED_STD_VECTOR(NodeModel_t)
        node_models_; // node models

        State_t state_desc_;     // state description
        Control_t control_desc_; // control description

        Scalar period_;    // time period
        size_t num_nodes_; // number of nodes
        size_t nc_;        // number of integration constraints

    }; // class SegmentModelBase

} // namespace galileo

#endif // __galileo_core_segment_model_base_hpp__
