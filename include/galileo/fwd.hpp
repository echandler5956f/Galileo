#ifndef __galileo_fwd_hpp__
#define __galileo_fwd_hpp__

// Forward declaration of the galileo namespace
namespace galileo
{
} // namespace galileo

#include <cassert>
#include <memory>
#include <type_traits>

#include "galileo/macros.hpp"

namespace galileo
{

    template <class T>
    struct traits
    {
    };

    namespace internal
    {

        struct Blank
        {
        }; // struct Blank

        // Base class for numerical classes.
        template <class Derived>
        struct CRTP
        {

            /** Return reference to this as derived object */
            inline Derived &derived() & noexcept
            {
                return *static_cast<Derived *>(this);
            }
            /** Return reference to this as derived object */
            inline const Derived &derived() const & noexcept
            {
                return *static_cast<const Derived *>(this);
            }
            /** Return reference to this as derived object, when this is rvalue */
            inline Derived &&derived() && noexcept
            {
                return std::move(*static_cast<Derived *>(this));
            }

        }; // struct CRTP

    } // namespace internal

} // namespace galileo

#include "galileo/eigen.hpp"

#endif // __galileo_fwd_hpp__
