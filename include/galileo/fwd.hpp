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
                return *static_cast<const Derived *>(this);
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

    // Helper alias that automatically corrects the storage order for 1xN or Nx1
    template <typename Scalar, int Rows, int Cols, int DefaultOptions>
    struct Matrix
    {
        static constexpr int OptionsFixed =
            (Rows == 1 && Cols != 1) ? (DefaultOptions | Eigen::RowMajor) : (Cols == 1 && Rows != 1) ? (DefaultOptions & ~Eigen::RowMajor)
                                                                                                        : DefaultOptions;

        using type = Eigen::Matrix<Scalar, Rows, Cols, OptionsFixed>;
    };

} // namespace galileo

#endif // __galileo_fwd_hpp__