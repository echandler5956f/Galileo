#ifndef __galileo_common_meta_dimensions_hpp__
#define __galileo_common_meta_dimensions_hpp__

#include <cassert>
#include <cmath>
#include <type_traits>

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
    template <int Value_ = detail::Dynamic>
    class Dimension
    {
    public:
        static constexpr int Value = Value_;
        static constexpr bool IsDynamic = (Value == detail::Dynamic);
        static constexpr bool IsFixed = !IsDynamic;

    private:
        int runtime_value_;

    public:
        constexpr Dimension()
            requires IsFixed
            : runtime_value_(Value)
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
            assert(runtime_val == Value &&
                   "Dimension value does not match fixed compile-time size");
        }

        // Copy and assignment
        Dimension(const Dimension &) = default;
        Dimension &operator=(const Dimension &) = default;

        // Value access
        constexpr int value() const { return runtime_value_; }
        constexpr operator int() const { return value(); }

        void set_value(int runtime_val)
        {
            if constexpr (IsDynamic)
                runtime_value_ = runtime_val;
            else
                assert(runtime_val == Value &&
                       "Dimension value does not match fixed compile-time size");
        }

        // Arithmetic operations
        template <int OtherValue>
        auto operator+(const Dimension<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::add<Value, OtherValue>::Value;
            return Dimension<result_compile_time>(value() + other.value());
        }

        template <int OtherValue>
        auto operator-(const Dimension<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::sub<Value, OtherValue>::Value;
            return Dimension<result_compile_time>(value() - other.value());
        }

        template <int OtherValue>
        auto operator*(const Dimension<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::mul<Value, OtherValue>::Value;
            return Dimension<result_compile_time>(value() * other.value());
        }

        template <int OtherValue>
        auto operator/(const Dimension<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::div<Value, OtherValue>::Value;
            return Dimension<result_compile_time>(value() / other.value());
        }

        // Scalar operations
        auto operator+(int scalar) const
        {
            constexpr int result_compile_time = IsFixed ? Value + scalar : detail::Dynamic;
            return Dimension<result_compile_time>(value() + scalar);
        }

        auto operator-(int scalar) const
        {
            constexpr int result_compile_time = IsFixed ? Value - scalar : detail::Dynamic;
            return Dimension<result_compile_time>(value() - scalar);
        }

        auto operator*(int scalar) const
        {
            constexpr int result_compile_time = IsFixed ? Value * scalar : detail::Dynamic;
            return Dimension<result_compile_time>(value() * scalar);
        }

        auto operator/(int scalar) const
        {
            constexpr int result_compile_time = IsFixed ? Value / scalar : detail::Dynamic;
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

    // Helper to extract compile-time value from raw integral or dimension types
    template <auto Val>
    struct extract_compile_time_value
    {
        static constexpr int Value = []()
        {
            if constexpr (std::is_integral_v<decltype(Val)>)
            {
                return Val;
            }
            else
            {
                static_assert(!Val.IsDynamic, "Cannot use dynamic dimension as template parameter");
                return Val.Value;
            }
        }();
    };

} // namespace galileo

#endif // __galileo_common_meta_dimensions_hpp__
