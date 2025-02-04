#ifndef __galileo_fwd_hpp__
#define __galileo_fwd_hpp__

// Forward declaration of the galileo namespace
namespace galileo
{
} // namespace galileo

#include <cassert>
#include <type_traits>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace galileo
{

    namespace internal
    {

        template <typename T>
        struct traits; // forward declaration

        // From Eigen:
        // Here we say once and for all that traits<const T> == traits<T>
        // When constness must affect traits, it has to be constness on template parameters on which T itself depends.
        // For example, traits<Map<const T> > != traits<Map<T> >, but
        //              traits<const Map<T> > == traits<Map<T> >
        template <typename T>
        struct traits<const T> : traits<T>
        {
        }; // specialization for const types

        struct Blank
        {
        }; // struct Blank

        // Base class for numerical classes.
        template <class Derived>
        struct CRTP
        {
            using NumScalar = typename traits<Derived>::NumScalar;
            using VarScalar = typename traits<Derived>::VarScalar;

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

        }; // struct CRTP

    } // namespace internal

} // namespace galileo

#endif // __galileo_fwd_hpp__