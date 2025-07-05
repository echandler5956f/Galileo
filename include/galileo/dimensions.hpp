#ifndef __galileo_dimensions_hpp__
#define __galileo_dimensions_hpp__

#include <cassert>
#include <cmath>

namespace galileo
{

    // Compile-time arithmetic helpers
    namespace detail
    {
        constexpr int Dynamic = -1; // Most sensibly, this should be the Eigen::Dynamic value (-1)

        // Addition with dynamic handling
        template <int A, int B>
        struct add
        {
            static constexpr int Value = (A == Dynamic || B == Dynamic) ? Dynamic : (A + B);
        };

        // Subtraction with dynamic handling
        template <int A, int B>
        struct sub
        {
            static constexpr int Value = (A == Dynamic || B == Dynamic) ? Dynamic : (A - B);
        };

        // Multiplication with dynamic handling
        template <int A, int B>
        struct mul
        {
            static constexpr int Value = (A == Dynamic || B == Dynamic) ? Dynamic : (A * B);
        };

        // Division with dynamic handling
        template <int A, int B>
        struct div
        {
            static constexpr int Value = (A == Dynamic || B == Dynamic || B == 0) ? Dynamic : (A / B);
        };

        // Maximum with dynamic handling
        template <int A, int B>
        struct max
        {
            static constexpr int Value = (A == Dynamic || B == Dynamic) ? Dynamic : ((A > B) ? A : B);
        };

        // Minimum with dynamic handling
        template <int A, int B>
        struct min
        {
            static constexpr int Value = (A == Dynamic || B == Dynamic) ? Dynamic : ((A < B) ? A : B);
        };

    } // namespace detail

    /**
     * @brief Unified dimension class that handles both compile-time and runtime dimensions
     *
     * This class encapsulates a dimension that may be either fixed at compile-time or
     * dynamic (known only at runtime). It provides seamless arithmetic operations and
     * automatic dispatch for Eigen operations.
     */
    template <int CompileTimeValue_ = detail::Dynamic>
    class Dimension
    {
    public:
        static constexpr int CompileTimeValue = CompileTimeValue_;
        static constexpr bool IsDynamic = (CompileTimeValue == detail::Dynamic);
        static constexpr bool IsFixed = !IsDynamic;

    private:
        int runtime_value_;

    public:
        constexpr Dimension()
            requires IsFixed
            : runtime_value_(CompileTimeValue)
        {
        }

        constexpr Dimension()
            requires IsDynamic
            : runtime_value_(0)
        {
        }

        constexpr Dimension(int runtime_val)
            requires IsDynamic
            : runtime_value_(runtime_val)
        {
        }

        constexpr explicit Dimension(int runtime_val)
            requires(!IsDynamic)
            : runtime_value_(runtime_val)
        {
            assert(runtime_val == CompileTimeValue &&
                   "Dimension value does not match fixed compile-time size");
        }

        // Copy and assignment
        Dimension(const Dimension &) = default;
        Dimension &operator=(const Dimension &) = default;

        // Conversion from fixed to dynamic dimension
        template <int OtherValue>
        Dimension(const Dimension<OtherValue> &other)
            : runtime_value_(other.value()) {}

        // Value access
        constexpr int value() const { return runtime_value_; }
        constexpr operator int() const { return value(); }

        // Compile-time value access (for template parameters)
        static constexpr int cvalue()
        {
            return CompileTimeValue;
        }

        // Arithmetic operations
        template <int OtherValue>
        auto operator+(const Dimension<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::add<CompileTimeValue, OtherValue>::Value;
            return Dimension<result_compile_time>(value() + other.value());
        }

        template <int OtherValue>
        auto operator-(const Dimension<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::sub<CompileTimeValue, OtherValue>::Value;
            return Dimension<result_compile_time>(value() - other.value());
        }

        template <int OtherValue>
        auto operator*(const Dimension<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::mul<CompileTimeValue, OtherValue>::Value;
            return Dimension<result_compile_time>(value() * other.value());
        }

        template <int OtherValue>
        auto operator/(const Dimension<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::div<CompileTimeValue, OtherValue>::Value;
            return Dimension<result_compile_time>(value() / other.value());
        }

        // Scalar operations
        auto operator+(int scalar) const
        {
            constexpr int result_compile_time = IsFixed ? CompileTimeValue + scalar : detail::Dynamic;
            return Dimension<result_compile_time>(value() + scalar);
        }

        auto operator-(int scalar) const
        {
            constexpr int result_compile_time = IsFixed ? CompileTimeValue - scalar : detail::Dynamic;
            return Dimension<result_compile_time>(value() - scalar);
        }

        auto operator*(int scalar) const
        {
            constexpr int result_compile_time = IsFixed ? CompileTimeValue * scalar : detail::Dynamic;
            return Dimension<result_compile_time>(value() * scalar);
        }

        auto operator/(int scalar) const
        {
            constexpr int result_compile_time = IsFixed ? CompileTimeValue / scalar : detail::Dynamic;
            return Dimension<result_compile_time>(value() / scalar);
        }

        // Comparison operations
        bool operator==(const Dimension &other) const { return value() == other.value(); }
        bool operator!=(const Dimension &other) const { return value() != other.value(); }
        bool operator<(const Dimension &other) const { return value() < other.value(); }
        bool operator<=(const Dimension &other) const { return value() <= other.value(); }
        bool operator>(const Dimension &other) const { return value() > other.value(); }
        bool operator>=(const Dimension &other) const { return value() >= other.value(); }
    }; // class Dimension

    // Type aliases for common cases
    using Dynamic = Dimension<detail::Dynamic>;
    template <int N>
    using Fixed = Dimension<N>;

    // Factory functions
    inline Dynamic dynamic(int runtime_value) { return Dynamic(runtime_value); }

    template <int N>
    inline Fixed<N> fixed() { return Fixed<N>(); }

    // Maximum and minimum operations
    template <int A, int B>
    auto max(const Dimension<A> &a, const Dimension<B> &b)
    {
        constexpr int result_compile_time = detail::max<A, B>::Value;
        return Dimension<result_compile_time>(std::max(a.value(), b.value()));
    }

    template <int A, int B>
    auto min(const Dimension<A> &a, const Dimension<B> &b)
    {
        constexpr int result_compile_time = detail::min<A, B>::Value;
        return Dimension<result_compile_time>(std::min(a.value(), b.value()));
    }

} // namespace galileo

#endif // __galileo_dimensions_hpp__
