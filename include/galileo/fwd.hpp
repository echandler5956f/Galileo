#pragma once

namespace galileo
{
}

#include <Eigen/Core>

namespace galileo
{

    // Common traits structure to fully define base classes for CRTP.
    template <class C>
    struct traits
    {
    };

    template <class Derived>
    struct NumericalBase
    {
        typedef typename traits<Derived>::Scalar Scalar;
    };

    // Type of the cast of a class C templated by Scalar and Options, to a new NewScalar type.
    // This class should be specialized for each types.
    template <typename NewScalar, class C>
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
}