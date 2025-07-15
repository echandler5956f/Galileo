#ifndef __galileo_common_meta_macros_hpp__
#define __galileo_common_meta_macros_hpp__

#define GALILEO_WORLD_VERSION 2
#define GALILEO_MAJOR_VERSION 0
#define GALILEO_MINOR_VERSION 00

#define GALILEO_VERSION_AT_LEAST(x, y, z) (GALILEO_WORLD_VERSION > x || (GALILEO_WORLD_VERSION >= x &&                                \
                                                                         (GALILEO_MAJOR_VERSION > y || (GALILEO_MAJOR_VERSION >= y && \
                                                                                                        GALILEO_MINOR_VERSION >= z))))

// Custom assertion macro that can be configured for testing
#ifndef GALILEO_ASSERT
#ifdef GALILEO_TESTING
#include <stdexcept>
#define GALILEO_ASSERT(condition, message)                          \
    do                                                              \
    {                                                               \
        if (!(condition))                                           \
        {                                                           \
            throw std::runtime_error("Assertion failed: " message); \
        }                                                           \
    } while (0)
#else
#include <cassert>
#define GALILEO_ASSERT(condition, message) assert((condition) && (message))
#endif
#endif

/**
 * Forward accessor macro.
 *
 * This macro is used to forward the accessor to the derived class.
 *
 * @param ReturnType The return type of the accessor.
 * @param accessor_name The name of the accessor.
 */
#define FORWARD_ACCESSOR(ReturnType, accessor_name)        \
    /* lvalue-qualified overload */                        \
    ReturnType &accessor_name()                            \
    {                                                      \
        return this->derived().accessor_name##_accessor(); \
    }                                                      \
    /* const-lvalue-qualified overload */                  \
    const ReturnType &accessor_name() const                \
    {                                                      \
        return this->derived().accessor_name##_accessor(); \
    }

/**
 * Default accessor macro.
 *
 * This macro is used to define a default accessor for a member variable. This accessor macro is used in the
 * derived class that actually stores the accessor_name variable.
 *
 * @param ReturnType The return type of the accessor.
 * @param accessor_name The name of the accessor (name of the variable).
 */
#define DEFAULT_ACCESSOR(ReturnType, accessor_name)    \
    ReturnType &accessor_name##_accessor()             \
    {                                                  \
        return accessor_name;                          \
    }                                                  \
    const ReturnType &accessor_name##_accessor() const \
    {                                                  \
        return accessor_name;                          \
    }

/**
 * Generic accessor macro for type-erased access to CRTP derived objects.
 *
 * This macro is used to define a generic accessor for a member variable (e.g., for a type-erased or Variant class).
 * See `include/galileo/core/costs/cost-generic.hpp`, `include/galileo/core/constraints/constraint-generic.hpp`,
 * `include/galileo/multibody/contacts/contact-generic.hpp`, and `include/galileo/predictive/phases/phase-generic.hpp`
 * for examples.
 *
 * @param ReturnType The return type of the accessor.
 * @param accessor_name The name of the accessor.
 */
#define GENERIC_ACCESSOR(ReturnType, accessor_name)    \
    ReturnType &accessor_name##_accessor()             \
    {                                                  \
        return accessor_name();                        \
    }                                                  \
    const ReturnType &accessor_name##_accessor() const \
    {                                                  \
        return accessor_name();                        \
    }

#endif // __galileo_common_meta_macros_hpp__
