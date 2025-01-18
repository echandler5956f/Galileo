#ifndef __galileo_core_segment_data_base_hpp__
#define __galileo_core_segment_data_base_hpp__

#include "galileo/core/segment/segment-base.hpp"
#include "galileo/core/segment/segment-model-base.hpp"

#include "galileo/utils/aligned-vector.hpp"

#define GALILEO_SEGMENT_DATA_TYPEDEF_GENERIC(Segment, TYPENAME)                                \
    GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, TYPENAME);                                  \
    typedef TYPENAME traits<Segment>::NodesDataVectorTypeConstRef NodesDataVectorTypeConstRef; \
    typedef TYPENAME traits<Segment>::NodesDataVectorTypeRef NodesDataVectorTypeRef;           \
    typedef TYPENAME traits<Segment>::StateTypeConstRef StateTypeConstRef;                     \
    typedef TYPENAME traits<Segment>::StateTypeRef StateTypeRef;                               \
    typedef TYPENAME traits<Segment>::ControlTypeConstRef ControlTypeConstRef;                 \
    typedef TYPENAME traits<Segment>::ControlTypeRef ControlTypeRef;                           \
    typedef TYPENAME traits<Segment>::CTypeConstRef CTypeConstRef;                             \
    typedef TYPENAME traits<Segment>::CTypeRef CTypeRef;                                       \
    typedef TYPENAME traits<Segment>::CkTypeConstRef CkTypeConstRef;                           \
    typedef TYPENAME traits<Segment>::CkTypeRef CkTypeRef;                                     \
    typedef TYPENAME traits<Segment>::CwTypeConstRef CwTypeConstRef;                           \
    typedef TYPENAME traits<Segment>::CwTypeRef CwTypeRef;

#ifdef __clang__

#define GALILEO_SEGMENT_DATA_TYPEDEF(Segment) \
    GALILEO_SEGMENT_DATA_TYPEDEF_GENERIC(Segment, GALILEO_EMPTY_ARG)
#define GALILEO_SEGMENT_DATA_TYPEDEF_TEMPLATE(Segment) \
    GALILEO_SEGMENT_DATA_TYPEDEF_GENERIC(Segment, typename)

#elif (__GNUC__ == 4) && (__GNUC_MINOR__ == 4) && (__GNUC_PATCHLEVEL__ == 2)

#define GALILEO_SEGMENT_DATA_TYPEDEF(Segment) \
    GALILEO_SEGMENT_DATA_TYPEDEF_GENERIC(Segment, GALILEO_EMPTY_ARG)
#define GALILEO_SEGMENT_DATA_TYPEDEF_TEMPLATE(Segment) \
    GALILEO_SEGMENT_DATA_TYPEDEF_GENERIC(Segment, typename)

#else

#define GALILEO_SEGMENT_DATA_TYPEDEF(Segment) GALILEO_SEGMENT_DATA_TYPEDEF_GENERIC(Segment, typename)
#define GALILEO_SEGMENT_DATA_TYPEDEF_TEMPLATE(Segment) \
    GALILEO_SEGMENT_DATA_TYPEDEF_GENERIC(Segment, typename)

#endif

#define GALILEO_SEGMENT_DATA_BASE_DEFAULT_ACCESSOR          \
    NodesDataVectorTypeConstRef nodes_data_accessor() const \
    {                                                       \
        return nodes_data;                                  \
    }                                                       \
    NodesDataVectorTypeRef nodes_data_accessor()            \
    {                                                       \
        return nodes_data;                                  \
    }                                                       \
    CTypeConstRef C_accessor() const                        \
    {                                                       \
        return C;                                           \
    }                                                       \
    CTypeRef C_accessor()                                   \
    {                                                       \
        return C;                                           \
    }                                                       \
    CkTypeConstRef Ck_accessor() const                      \
    {                                                       \
        return Ck;                                          \
    }                                                       \
    CkTypeRef Ck_accessor()                                 \
    {                                                       \
        return Ck;                                          \
    }                                                       \
    CwTypeConstRef Cw_accessor() const                      \
    {                                                       \
        return Cw;                                          \
    }                                                       \
    CwTypeRef Cw_accessor()                                 \
    {                                                       \
        return Cw;                                          \
    }

#define GALILEO_SEGMENT_DATA_BASE_ACCESSOR_DEFAULT_RETURN_TYPE                          \
    typedef const GALILEO_ALIGNED_STD_VECTOR(NodeData_t) & NodesDataVectorTypeConstRef; \
    typedef GALILEO_ALIGNED_STD_VECTOR(NodeData_t) & NodesDataVectorTypeRef;            \
    typedef const State_t &StateTypeConstRef;                                           \
    typedef State_t &StateTypeRef;                                                      \
    typedef const Control_t &ControlTypeConstRef;                                       \
    typedef Control_t &ControlTypeRef;                                                  \
    typedef const C_t &CTypeConstRef;                                                   \
    typedef C_t &CTypeRef;                                                              \
    typedef const Eigen::Matrix<Scalar, Nc, Nk> &CkTypeConstRef;                        \
    typedef Eigen::Matrix<Scalar, Nc, Nk> &CkTypeRef;                                   \
    typedef const Eigen::Matrix<Scalar, Nc, Nw> &CwTypeConstRef;                        \
    typedef Eigen::Matrix<Scalar, Nc, Nw> &CwTypeRef;

namespace galileo
{

    template <typename Derived>
    struct SegmentDataBase : CRTP<Derived>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using SegmentDerived = typename traits<Derived>::SegmentDerived;
        GALILEO_SEGMENT_DATA_TYPEDEF_TEMPLATE(SegmentDerived);

        NodesDataVectorTypeConstRef nodes_data() const
        {
            return derived().nodes_data_accessor();
        }
        NodesDataVectorTypeRef nodes_data()
        {
            return derived().nodes_data_accessor();
        }

        CTypeConstRef C() const
        {
            return derived().C_accessor();
        }
        CTypeRef C()
        {
            return derived().C_accessor();
        }

        CkTypeConstRef Ck() const
        {
            return derived().Ck_accessor();
        }
        CkTypeRef Ck()
        {
            return derived().Ck_accessor();
        }

        CwTypeConstRef Cw() const
        {
            return derived().Cw_accessor();
        }
        CwTypeRef Cw()
        {
            return derived().Cw_accessor();
        }

        // template <typename OtherDerived>
        // bool operator==(const SegmentDataBase<OtherDerived> &other) const
        // {
        //     return derived().isEqual(other.derived());
        // }

        // // Default operator== implementation
        // bool isEqual(const SegmentDataBase<Derived> & other) const
        // {
        // return internal::comparison_eq(joint_q(), other.joint_q())
        //         && internal::comparison_eq(joint_v(), other.joint_v())
        //         && internal::comparison_eq(S(), other.S()) && internal::comparison_eq(M(), other.M())
        //         && internal::comparison_eq(v(), other.v()) && internal::comparison_eq(c(), other.c())
        //         && internal::comparison_eq(U(), other.U())
        //         && internal::comparison_eq(Dinv(), other.Dinv())
        //         && internal::comparison_eq(UDinv(), other.UDinv());
        // }

        // // Default operator== implementation
        // template <typename OtherDerived>
        // bool isEqual(const SegmentDataBase<OtherDerived> & /*other*/) const
        // {
        //     return false;
        // }

        // bool operator!=(const SegmentDataBase<Derived> &other) const
        // {
        //     return derived().isNotEqual(other.derived());
        // }

        // // Default operator!= implementation
        // bool isNotEqual(const SegmentDataBase<Derived> &other) const
        // {
        //     return !(internal::comparison_eq(derived(), other.derived()));
        // }

    protected:
        // Default constructor: protected.
        inline SegmentDataBase()
        {
        }

    }; // struct SegmentDataBase

} // namespace galileo

#endif // __galileo_core_segment_data_base_hpp__