#ifndef __galileo_fwd_hpp__
#define __galileo_fwd_hpp__

namespace galileo
{
} // namespace galileo

#include <cassert>
#include <type_traits>

#include <Eigen/Core>
#include "galileo/utils/cast.hpp"
#include "galileo/utils/macro.hpp"


namespace galileo
{
    /**
     * @brief The kalmanif traits class.
     *
     * @tparam T The type for which to specialize the traits
     * @tparam Enable a SFINAE argument
     */
    template <typename T, class Enable = void>
    struct traits;

    /// @note the following is from the Eigen library
    /// here we say once and for all that traits<const T> == traits<T>
    ///
    /// When constness must affect traits, it has to be constness on
    /// template parameters on which T itself depends.
    /// For example, traits<Map<const T> > != traits<Map<T> >, but
    ///              traits<const Map<T> > == traits<Map<T> >
    template <typename T, class Enable>
    struct traits<const T, Enable> : traits<T, Enable>
    {
    };

    // Blank type
    struct Blank
    {
    };

    // Base class for numerical classes.
    template <class Derived>
    struct CRTP
    {
        using Scalar = typename traits<Derived>::Scalar;

    protected:
        /** Return reference to this as derived object */
        inline Derived &derived() & noexcept
        {
            return *static_cast<Derived *>(this);
        }
        /** Return reference to this as derived object */
        inline const Derived &derived() const & noexcept
        {
            return *static_cast<Derived const *>(this);
        }
        /** Return reference to this as derived object, when this is rvalue */
        inline Derived &&derived() && noexcept
        {
            return std::move(*static_cast<Derived *>(this));
        }
    };

    // Type of the cast of a class T templated by Scalar and Options, to a new NewScalar type.
    // This class should be specialized for each types.
    template <typename NewScalar, class T>
    struct CastType;

    // Cast scalar type from type FROM to type TO.
    template <typename To, typename From>
    struct ScalarCast
    {
        static To cast(const From &value)
        {
            return static_cast<To>(value);
        }
    };

    template <typename To, typename From>
    To scalar_cast(const From &value)
    {
        return ScalarCast<To, From>::cast(value);
    }

    enum AssignmentOp
    {
        setto,
        addto,
        rmfrom
    };

    inline bool is_a_AssignmentOp(AssignmentOp op)
    {
        return (op == setto || op == addto || op == rmfrom);
    }

    struct ReturnTypeNotDefined;

} // namespace galileo

#include "galileo/context.hpp"

#endif // __galileo_fwd_hpp__