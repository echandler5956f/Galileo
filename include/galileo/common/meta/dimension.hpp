#ifndef __galileo_common_meta_dimension_hpp__
#define __galileo_common_meta_dimension_hpp__

#include "galileo/common/meta/macros.hpp"

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
            static constexpr int Value = (A == Dynamic || B == Dynamic || A < B) ? Dynamic : (A - B);
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
    class DimensionTpl
    {
    public:
        static constexpr int Value = Value_;
        static constexpr bool IsDynamic = (Value == detail::Dynamic);
        static constexpr bool IsFixed = !IsDynamic;

    private:
        int runtime_value_;

    public:
        constexpr DimensionTpl()
            requires IsFixed
            : runtime_value_(Value)
        {
            static_assert(Value >= 0, "Compile-time dimension values must be non-negative");
        }

        constexpr DimensionTpl()
            requires IsDynamic
            : runtime_value_(0)
        {
        }

        constexpr DimensionTpl(int runtime_val)
            requires IsDynamic
            : runtime_value_(runtime_val)
        {
            GALILEO_ASSERT(runtime_val >= 0, "DimensionTpl: Runtime dimension values must be non-negative");
        }

        constexpr explicit DimensionTpl(int runtime_val)
            requires(!IsDynamic)
            : runtime_value_(runtime_val)
        {
            GALILEO_ASSERT(runtime_val == Value, "DimensionTpl: Dimension value does not match fixed compile-time size");
        }

        // Copy and assignment
        DimensionTpl(const DimensionTpl &) = default;
        DimensionTpl &operator=(const DimensionTpl &) = default;

        // Value access
        constexpr int value() const { return runtime_value_; }
        constexpr operator int() const { return value(); }

        void set_value(int runtime_val)
        {
            if constexpr (IsDynamic)
            {
                GALILEO_ASSERT(runtime_val >= 0, "DimensionTpl: Runtime dimension values must be non-negative");
                runtime_value_ = runtime_val;
            }
            else
            {
                GALILEO_ASSERT(runtime_val == Value, "DimensionTpl: Dimension value does not match fixed compile-time size");
                GALILEO_ASSERT(runtime_val >= 0, "DimensionTpl: Runtime dimension values must be non-negative");
            }
        }

        // Arithmetic operations
        template <int OtherValue>
        auto operator+(const DimensionTpl<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::add<Value, OtherValue>::Value;
            int runtime_result = value() + other.value();
            GALILEO_ASSERT(runtime_result >= 0, "DimensionTpl: Dimension addition resulted in negative runtime value");
            return DimensionTpl<result_compile_time>(runtime_result);
        }

        template <int OtherValue>
        auto operator-(const DimensionTpl<OtherValue> &other) const
            requires(IsDynamic || OtherValue == detail::Dynamic || Value >= OtherValue)
        {
            constexpr int result_compile_time = detail::sub<Value, OtherValue>::Value;
            int runtime_result = value() - other.value();
            GALILEO_ASSERT(runtime_result >= 0, "DimensionTpl: Dimension subtraction resulted in negative runtime value");
            return DimensionTpl<result_compile_time>(runtime_result);
        }

        template <int OtherValue>
        auto operator*(const DimensionTpl<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::mul<Value, OtherValue>::Value;
            int runtime_result = value() * other.value();
            GALILEO_ASSERT(runtime_result >= 0, "DimensionTpl: Dimension multiplication resulted in negative runtime value");
            return DimensionTpl<result_compile_time>(runtime_result);
        }

        template <int OtherValue>
        auto operator/(const DimensionTpl<OtherValue> &other) const
        {
            constexpr int result_compile_time = detail::div<Value, OtherValue>::Value;
            int runtime_result = value() / other.value();
            GALILEO_ASSERT(runtime_result >= 0, "DimensionTpl: Dimension division resulted in negative runtime value");
            return DimensionTpl<result_compile_time>(runtime_result);
        }

        // Scalar operations
        auto operator+(int scalar) const
        {
            int runtime_result = value() + scalar;
            GALILEO_ASSERT(runtime_result >= 0, "DimensionTpl: Dimension scalar addition resulted in negative runtime value");
            return DimensionTpl<detail::Dynamic>(runtime_result);
        }

        auto operator-(int scalar) const
        {
            int runtime_result = value() - scalar;
            GALILEO_ASSERT(runtime_result >= 0, "DimensionTpl: Dimension scalar subtraction resulted in negative runtime value");
            return DimensionTpl<detail::Dynamic>(runtime_result);
        }

        auto operator*(int scalar) const
        {
            GALILEO_ASSERT(scalar >= 0, "DimensionTpl: Dimension scalar multiplication resulted in negative runtime value");
            return DimensionTpl<detail::Dynamic>(value() * scalar);
        }

        auto operator/(int scalar) const
        {
            GALILEO_ASSERT(scalar >= 0, "DimensionTpl: Dimension scalar division resulted in negative runtime value");
            return DimensionTpl<detail::Dynamic>(value() / scalar);
        }

        // Comparison operations
        bool operator==(const DimensionTpl &other) const { return value() == other.value(); }
        bool operator!=(const DimensionTpl &other) const { return value() != other.value(); }
        bool operator<(const DimensionTpl &other) const { return value() < other.value(); }
        bool operator<=(const DimensionTpl &other) const { return value() <= other.value(); }
        bool operator>(const DimensionTpl &other) const { return value() > other.value(); }
        bool operator>=(const DimensionTpl &other) const { return value() >= other.value(); }
    }; // class Dimension

    // Maximum and minimum operations
    template <int A, int B>
    auto max(const DimensionTpl<A> &a, const DimensionTpl<B> &b)
    {
        constexpr int result_compile_time = detail::max<A, B>::Value;
        return DimensionTpl<result_compile_time>(std::max(a.value(), b.value()));
    }

    template <int A, int B>
    auto min(const DimensionTpl<A> &a, const DimensionTpl<B> &b)
    {
        constexpr int result_compile_time = detail::min<A, B>::Value;
        return DimensionTpl<result_compile_time>(std::min(a.value(), b.value()));
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

#endif // __galileo_common_meta_dimension_hpp__
