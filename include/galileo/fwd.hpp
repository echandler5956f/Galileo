#ifndef __galileo_fwd_hpp__
#define __galileo_fwd_hpp__

// Forward declaration of the galileo namespace
namespace galileo
{
} // namespace galileo

#include <cassert>
#include <type_traits>
#include <memory>

#include <galileo/macros.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

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
                return *static_cast<Derived const *>(this);
            }
            /** Return reference to this as derived object, when this is rvalue */
            inline Derived &&derived() && noexcept
            {
                return std::move(*static_cast<Derived *>(this));
            }

        }; // struct CRTP

    } // namespace internal

    namespace compile_time
    {

        // We reserve -1 to represent a dynamic size.
        constexpr int dynamic_size = -1;

        // A strongly typed helper that encapsulates a compile-time size as a non-type template parameter.
        // Each operator overload is marked consteval to force compile-time evaluation.
        template <int N>
        struct Size
        {
            static constexpr int value = N;

            // Addition operator overload.
            template <int M>
            consteval auto operator+(Size<M>) const
            {
                return Size<(N == dynamic_size || M == dynamic_size ? dynamic_size : N + M)>{};
            }

            // Subtraction operator overload.
            template <int M>
            consteval auto operator-(Size<M>) const
            {
                return Size<(N == dynamic_size || M == dynamic_size ? dynamic_size : N - M)>{};
            }

            // Multiplication operator overload.
            template <int M>
            consteval auto operator*(Size<M>) const
            {
                return Size<(N == dynamic_size || M == dynamic_size ? dynamic_size : N * M)>{};
            }

            // Division operator overload.
            // The static_assert ensures division by zero is caught at compile time.
            template <int M>
            consteval auto operator/(Size<M>) const
            {
                static_assert(M != 0, "Division by zero is not allowed.");
                return Size<(N == dynamic_size || M == dynamic_size ? dynamic_size : N / M)>{};
            }
        };

    } // namespace compile_time

    enum AssignmentOp
    {
        setto,
        addto,
        rmfrom
    }; // enum AssignmentOp

    inline bool is_a_AssignmentOp(AssignmentOp op)
    {
        return (op == setto || op == addto || op == rmfrom);
    }

    enum Jcomponent
    {
        both = 0,
        first = 1,
        second = 2
    }; // enum Jcomponent

    inline bool is_a_Jcomponent(Jcomponent firstsecond)
    {
        return (firstsecond == first || firstsecond == second || firstsecond == both);
    }

    // Use an Eigen::Map if the matrix is shared.
    template <typename MatrixType, bool Share>
    using SelectMatrix = std::conditional_t<Share, Eigen::Map<MatrixType>, MatrixType>;

    //-----------------------------------------------------------------
    // Helper: Compute the index of a type T in a tuple.
    template <typename T, typename Tuple>
    struct tuple_index;

    template <typename T, typename... Ts>
    struct tuple_index<T, std::tuple<T, Ts...>>
    {
        static constexpr std::size_t value = 0;
    }; // struct tuple_index

    template <typename T, typename U, typename... Ts>
    struct tuple_index<T, std::tuple<U, Ts...>>
    {
        static constexpr std::size_t value = 1 + tuple_index<T, std::tuple<Ts...>>::value;
    }; // struct tuple_index

} // namespace galileo

#endif // __galileo_fwd_hpp__