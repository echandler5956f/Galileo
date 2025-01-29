#ifndef __galileo_core_phase_data_base_hpp__
#define __galileo_core_phase_data_base_hpp__

#include "galileo/core/phase/phase-base.hpp"
#include "galileo/core/phase/phase-model-base.hpp"

#include "galileo/utils/aligned-vector.hpp"

#define GALILEO_PHASE_DATA_TYPEDEF_GENERIC(Phase, TYPENAME)                                \
    GALILEO_PHASE_MODEL_TYPEDEF_GENERIC(Phase, TYPENAME);                                  \
    typedef TYPENAME traits<Phase>::NodesDataVectorTypeConstRef NodesDataVectorTypeConstRef; \
    typedef TYPENAME traits<Phase>::NodesDataVectorTypeRef NodesDataVectorTypeRef;           \
    typedef TYPENAME traits<Phase>::StateTypeConstRef StateTypeConstRef;                     \
    typedef TYPENAME traits<Phase>::StateTypeRef StateTypeRef;                               \
    typedef TYPENAME traits<Phase>::ControlTypeConstRef ControlTypeConstRef;                 \
    typedef TYPENAME traits<Phase>::ControlTypeRef ControlTypeRef;                           \
    typedef TYPENAME traits<Phase>::CTypeConstRef CTypeConstRef;                             \
    typedef TYPENAME traits<Phase>::CTypeRef CTypeRef;                                       \
    typedef TYPENAME traits<Phase>::CkTypeConstRef CkTypeConstRef;                           \
    typedef TYPENAME traits<Phase>::CkTypeRef CkTypeRef;                                     \
    typedef TYPENAME traits<Phase>::CwTypeConstRef CwTypeConstRef;                           \
    typedef TYPENAME traits<Phase>::CwTypeRef CwTypeRef;

#define GALILEO_PHASE_DATA_TYPEDEF_TEMPLATE(Phase) \
    GALILEO_PHASE_DATA_TYPEDEF_GENERIC(Phase, typename)

#define GALILEO_PHASE_DATA_BASE_DEFAULT_ACCESSOR          \
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

#define GALILEO_PHASE_DATA_BASE_ACCESSOR_DEFAULT_RETURN_TYPE                          \
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
    struct PhaseDataBase : CRTP<Derived>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PhaseDerived = typename traits<Derived>::PhaseDerived;
        GALILEO_PHASE_DATA_TYPEDEF_TEMPLATE(PhaseDerived);

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
        // bool operator==(const PhaseDataBase<OtherDerived> &other) const
        // {
        //     return derived().isEqual(other.derived());
        // }

        // // Default operator== implementation
        // bool isEqual(const PhaseDataBase<Derived> & other) const
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
        // bool isEqual(const PhaseDataBase<OtherDerived> & /*other*/) const
        // {
        //     return false;
        // }

        // bool operator!=(const PhaseDataBase<Derived> &other) const
        // {
        //     return derived().isNotEqual(other.derived());
        // }

        // // Default operator!= implementation
        // bool isNotEqual(const PhaseDataBase<Derived> &other) const
        // {
        //     return !(internal::comparison_eq(derived(), other.derived()));
        // }

    protected:
        // Default constructor: protected.
        inline PhaseDataBase()
        {
        }

    }; // struct PhaseDataBase

} // namespace galileo

#endif // __galileo_core_phase_data_base_hpp__