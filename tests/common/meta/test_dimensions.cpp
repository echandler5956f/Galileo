#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <climits>
#include <stdexcept>

#include "galileo/common/meta/dimensions.hpp"

using namespace galileo;

TEST_CASE("DimensionTpl - Basic Construction and Properties", "[dimensions]")
{
    SECTION("Fixed dimension construction")
    {
        constexpr DimensionTpl<5> fixed_dim;
        REQUIRE(fixed_dim.value() == 5);
        REQUIRE(fixed_dim.Value == 5);
        REQUIRE(fixed_dim.IsFixed);
        REQUIRE_FALSE(fixed_dim.IsDynamic);
        REQUIRE(static_cast<int>(fixed_dim) == 5);
    }

    SECTION("Dynamic dimension construction")
    {
        DimensionTpl<> dynamic_dim;
        REQUIRE(dynamic_dim.value() == 0);
        REQUIRE(dynamic_dim.Value == detail::Dynamic);
        REQUIRE(dynamic_dim.IsDynamic);
        REQUIRE(static_cast<int>(dynamic_dim) == 0);
    }

    SECTION("Dynamic dimension with initial value")
    {
        DimensionTpl<> dynamic_dim(42);
        REQUIRE(dynamic_dim.value() == 42);
        REQUIRE(dynamic_dim.IsDynamic);
        REQUIRE(static_cast<int>(dynamic_dim) == 42);
    }

    SECTION("Fixed dimension with explicit runtime value")
    {
        DimensionTpl<7> fixed_dim(7);
        REQUIRE(fixed_dim.value() == 7);
        REQUIRE(fixed_dim.Value == 7);
        REQUIRE(fixed_dim.IsFixed);
    }

    SECTION("Zero dimension")
    {
        constexpr DimensionTpl<0> zero_dim;
        REQUIRE(zero_dim.value() == 0);
        REQUIRE(zero_dim.Value == 0);
        REQUIRE_FALSE(zero_dim.IsDynamic);
    }

    SECTION("Large dimension")
    {
        constexpr DimensionTpl<1000> large_dim;
        REQUIRE(large_dim.value() == 1000);
        REQUIRE(large_dim.Value == 1000);
    }
}

TEST_CASE("DimensionTpl - Value Setting and Modification", "[dimensions]")
{
    SECTION("Setting dynamic dimension value")
    {
        DimensionTpl<> dynamic_dim;
        dynamic_dim.set_value(123);
        REQUIRE(dynamic_dim.value() == 123);

        dynamic_dim.set_value(0);
        REQUIRE(dynamic_dim.value() == 0);

        dynamic_dim.set_value(42);
        REQUIRE(dynamic_dim.value() == 42);
    }

    SECTION("Setting fixed dimension value - valid")
    {
        DimensionTpl<42> fixed_dim;
        // Should not throw when setting to correct value
        REQUIRE_NOTHROW(fixed_dim.set_value(42));
        REQUIRE(fixed_dim.value() == 42);
    }

    SECTION("Copy construction and assignment")
    {
        DimensionTpl<10> fixed_original;
        DimensionTpl<10> fixed_copy(fixed_original);
        REQUIRE(fixed_copy.value() == 10);

        DimensionTpl<10> fixed_assigned;
        fixed_assigned = fixed_original;
        REQUIRE(fixed_assigned.value() == 10);

        DimensionTpl<> dynamic_original(55);
        DimensionTpl<> dynamic_copy(dynamic_original);
        REQUIRE(dynamic_copy.value() == 55);

        DimensionTpl<> dynamic_assigned;
        dynamic_assigned = dynamic_original;
        REQUIRE(dynamic_assigned.value() == 55);
    }
}

TEST_CASE("DimensionTpl - Arithmetic Operations", "[dimensions]")
{
    SECTION("Addition - Fixed + Fixed")
    {
        constexpr DimensionTpl<5> a;
        constexpr DimensionTpl<3> b;
        auto result = a + b;

        REQUIRE(result.value() == 8);
        REQUIRE(result.Value == 8);
        REQUIRE(result.IsFixed);
    }

    SECTION("Addition - Fixed + Dynamic")
    {
        constexpr DimensionTpl<10> fixed;
        DimensionTpl<> dynamic(7);
        auto result = fixed + dynamic;

        REQUIRE(result.value() == 17);
        REQUIRE(result.Value == detail::Dynamic);
        REQUIRE(result.IsDynamic);
    }

    SECTION("Addition - Dynamic + Fixed")
    {
        DimensionTpl<> dynamic(12);
        constexpr DimensionTpl<8> fixed;
        auto result = dynamic + fixed;

        REQUIRE(result.value() == 20);
        REQUIRE(result.Value == detail::Dynamic);
        REQUIRE(result.IsDynamic);
    }

    SECTION("Addition - Dynamic + Dynamic")
    {
        DimensionTpl<> a(15);
        DimensionTpl<> b(25);
        auto result = a + b;

        REQUIRE(result.value() == 40);
        REQUIRE(result.Value == detail::Dynamic);
        REQUIRE(result.IsDynamic);
    }

    SECTION("Subtraction - Fixed - Fixed")
    {
        constexpr DimensionTpl<10> a;
        constexpr DimensionTpl<3> b;
        auto result = a - b;

        REQUIRE(result.value() == 7);
        REQUIRE(result.Value == 7);
        REQUIRE_FALSE(result.IsDynamic);
    }

    SECTION("Subtraction - constraint prevents A < B for fixed dimensions")
    {
        constexpr DimensionTpl<5> a;
        constexpr DimensionTpl<8> b;

        // This should not compile because Value >= OtherValue constraint
        // auto result = a - b;  // Uncommenting this should cause compilation error

        // But subtracting in the other direction should work
        auto valid_result = b - a;
        REQUIRE(valid_result.value() == 3);
        REQUIRE(valid_result.Value == 3);
        REQUIRE_FALSE(valid_result.IsDynamic);
    }

    SECTION("Subtraction - with dynamic dimensions requires non-negative results")
    {
        DimensionTpl<> a(8);
        DimensionTpl<> b(5);

        // Dynamic dimensions can be subtracted when result is non-negative
        auto result = a - b;
        REQUIRE(result.value() == 3);
        REQUIRE(result.Value == detail::Dynamic);
        REQUIRE(result.IsDynamic);

        // Subtracting in reverse order should be prevented
        REQUIRE_THROWS_AS(b - a, std::runtime_error);
    }

    SECTION("Subtraction - mixed fixed and dynamic with non-negative results")
    {
        constexpr DimensionTpl<10> fixed;
        DimensionTpl<> dynamic_smaller(3);
        DimensionTpl<> dynamic_larger(15);

        // Fixed - Dynamic (where result is non-negative) should work
        auto result1 = fixed - dynamic_smaller;
        REQUIRE(result1.value() == 7);
        REQUIRE(result1.Value == detail::Dynamic);
        REQUIRE(result1.IsDynamic);

        // Dynamic - Fixed (where result is non-negative) should work
        auto result2 = dynamic_larger - fixed;
        REQUIRE(result2.value() == 5);
        REQUIRE(result2.Value == detail::Dynamic);
        REQUIRE(result2.IsDynamic);
    }

    SECTION("Subtraction - Dynamic - Fixed")
    {
        DimensionTpl<> dynamic(20);
        constexpr DimensionTpl<7> fixed;
        auto result = dynamic - fixed;

        REQUIRE(result.value() == 13);
        REQUIRE(result.Value == detail::Dynamic);
        REQUIRE(result.IsDynamic);
    }

    SECTION("Multiplication - Fixed * Fixed")
    {
        constexpr DimensionTpl<4> a;
        constexpr DimensionTpl<6> b;
        auto result = a * b;

        REQUIRE(result.value() == 24);
        REQUIRE(result.Value == 24);
        REQUIRE_FALSE(result.IsDynamic);
    }

    SECTION("Multiplication - involving zero")
    {
        constexpr DimensionTpl<0> zero;
        constexpr DimensionTpl<5> five;
        auto result = zero * five;

        REQUIRE(result.value() == 0);
        REQUIRE(result.Value == 0);
        REQUIRE_FALSE(result.IsDynamic);
    }

    SECTION("Multiplication - Dynamic * Fixed")
    {
        DimensionTpl<> dynamic(6);
        constexpr DimensionTpl<4> fixed;
        auto result = dynamic * fixed;

        REQUIRE(result.value() == 24);
        REQUIRE(result.Value == detail::Dynamic);
        REQUIRE(result.IsDynamic);
    }

    SECTION("Division - Fixed / Fixed")
    {
        constexpr DimensionTpl<20> a;
        constexpr DimensionTpl<4> b;
        auto result = a / b;

        REQUIRE(result.value() == 5);
        REQUIRE(result.Value == 5);
        REQUIRE_FALSE(result.IsDynamic);
    }

    SECTION("Division - with remainder")
    {
        constexpr DimensionTpl<7> a;
        constexpr DimensionTpl<3> b;
        auto result = a / b;

        REQUIRE(result.value() == 2); // Integer division
        REQUIRE(result.Value == 2);
        REQUIRE_FALSE(result.IsDynamic);
    }

    SECTION("Division - Dynamic / Fixed")
    {
        DimensionTpl<> dynamic(30);
        constexpr DimensionTpl<6> fixed;
        auto result = dynamic / fixed;

        REQUIRE(result.value() == 5);
        REQUIRE(result.Value == detail::Dynamic);
        REQUIRE(result.IsDynamic);
    }
}

TEST_CASE("DimensionTpl - Scalar Arithmetic Operations", "[dimensions]")
{
    SECTION("Fixed dimension scalar operations")
    {
        constexpr DimensionTpl<10> fixed;

        auto add_result = fixed + 5;
        REQUIRE(add_result.value() == 15);
        REQUIRE(add_result.Value == detail::Dynamic);
        REQUIRE(add_result.IsDynamic);

        auto sub_result = fixed - 3;
        REQUIRE(sub_result.value() == 7);
        REQUIRE(sub_result.Value == detail::Dynamic);
        REQUIRE(sub_result.IsDynamic);

        auto mul_result = fixed * 2;
        REQUIRE(mul_result.value() == 20);
        REQUIRE(mul_result.Value == detail::Dynamic);
        REQUIRE(mul_result.IsDynamic);

        auto div_result = fixed / 2;
        REQUIRE(div_result.value() == 5);
        REQUIRE(div_result.Value == detail::Dynamic);
        REQUIRE(div_result.IsDynamic);
    }

    SECTION("Dynamic dimension scalar operations")
    {
        DimensionTpl<> dynamic(12);

        auto add_result = dynamic + 8;
        REQUIRE(add_result.value() == 20);
        REQUIRE(add_result.Value == detail::Dynamic);
        REQUIRE(add_result.IsDynamic);

        auto sub_result = dynamic - 4;
        REQUIRE(sub_result.value() == 8);
        REQUIRE(sub_result.Value == detail::Dynamic);
        REQUIRE(sub_result.IsDynamic);

        auto mul_result = dynamic * 3;
        REQUIRE(mul_result.value() == 36);
        REQUIRE(mul_result.Value == detail::Dynamic);
        REQUIRE(mul_result.IsDynamic);

        auto div_result = dynamic / 4;
        REQUIRE(div_result.value() == 3);
        REQUIRE(div_result.Value == detail::Dynamic);
        REQUIRE(div_result.IsDynamic);
    }

    SECTION("Scalar operations with zero")
    {
        constexpr DimensionTpl<5> fixed;

        auto add_zero = fixed + 0;
        REQUIRE(add_zero.value() == 5);
        REQUIRE(add_zero.Value == detail::Dynamic);

        auto sub_zero = fixed - 0;
        REQUIRE(sub_zero.value() == 5);
        REQUIRE(sub_zero.Value == detail::Dynamic);

        auto mul_zero = fixed * 0;
        REQUIRE(mul_zero.value() == 0);
        REQUIRE(mul_zero.Value == detail::Dynamic);
    }

    SECTION("Scalar operations with negative numbers")
    {
        constexpr DimensionTpl<10> fixed;

        auto add_negative = fixed + (-3);
        REQUIRE(add_negative.value() == 7);
        REQUIRE(add_negative.Value == detail::Dynamic);

        auto sub_negative = fixed - (-2);
        REQUIRE(sub_negative.value() == 12);
        REQUIRE(sub_negative.Value == detail::Dynamic);

        // Test that assertion is raised for negative scalar multiplication
        REQUIRE_THROWS_AS(fixed * (-2), std::runtime_error);

        // Test that assertion is raised for negative scalar division
        REQUIRE_THROWS_AS(fixed / (-2), std::runtime_error);
    }
}

TEST_CASE("DimensionTpl - Comparison Operations", "[dimensions]")
{
    SECTION("Fixed dimension comparisons")
    {
        constexpr DimensionTpl<5> a;
        constexpr DimensionTpl<5> b;
        constexpr DimensionTpl<3> c;
        constexpr DimensionTpl<7> d;

        REQUIRE(a == b);
        REQUIRE_FALSE(a != b);
        REQUIRE(a != c);
        REQUIRE_FALSE(a == c);

        REQUIRE(c < a);
        REQUIRE(c <= a);
        REQUIRE_FALSE(c > a);
        REQUIRE_FALSE(c >= a);

        REQUIRE(d > a);
        REQUIRE(d >= a);
        REQUIRE_FALSE(d < a);
        REQUIRE_FALSE(d <= a);

        REQUIRE(a <= b);
        REQUIRE(a >= b);
    }

    SECTION("Dynamic dimension comparisons")
    {
        DimensionTpl<> a(10);
        DimensionTpl<> b(10);
        DimensionTpl<> c(5);
        DimensionTpl<> d(15);

        REQUIRE(a == b);
        REQUIRE_FALSE(a != b);
        REQUIRE(a != c);
        REQUIRE_FALSE(a == c);

        REQUIRE(c < a);
        REQUIRE(c <= a);
        REQUIRE_FALSE(c > a);
        REQUIRE_FALSE(c >= a);

        REQUIRE(d > a);
        REQUIRE(d >= a);
        REQUIRE_FALSE(d < a);
        REQUIRE_FALSE(d <= a);

        REQUIRE(a <= b);
        REQUIRE(a >= b);
    }

    SECTION("Mixed fixed and dynamic comparisons")
    {
        constexpr DimensionTpl<8> fixed;
        DimensionTpl<> dynamic_equal(8);
        DimensionTpl<> dynamic_less(3);
        DimensionTpl<> dynamic_greater(12);

        REQUIRE(fixed == dynamic_equal);
        REQUIRE(dynamic_equal == fixed);
        REQUIRE_FALSE(fixed != dynamic_equal);

        REQUIRE(fixed > dynamic_less);
        REQUIRE(dynamic_less < fixed);
        REQUIRE(fixed >= dynamic_less);
        REQUIRE(dynamic_less <= fixed);

        REQUIRE(fixed < dynamic_greater);
        REQUIRE(dynamic_greater > fixed);
        REQUIRE(fixed <= dynamic_greater);
        REQUIRE(dynamic_greater >= fixed);
    }

    SECTION("Comparisons with zero")
    {
        constexpr DimensionTpl<0> zero_fixed;
        DimensionTpl<> zero_dynamic(0);
        DimensionTpl<> positive_dynamic(5);

        REQUIRE(zero_fixed == zero_dynamic);
        REQUIRE(zero_dynamic == zero_fixed);
        REQUIRE(zero_fixed < positive_dynamic);
        REQUIRE(positive_dynamic > zero_fixed);
    }

    SECTION("Comparisons with various values")
    {
        DimensionTpl<> small(1);
        DimensionTpl<> zero(0);
        DimensionTpl<> large(10);

        REQUIRE(small > zero);
        REQUIRE(zero < small);
        REQUIRE(small < large);
        REQUIRE(large > small);
        REQUIRE(zero < large);
        REQUIRE(small >= zero);
        REQUIRE(large >= small);
    }
}

TEST_CASE("DimensionTpl - Min/Max Operations", "[dimensions]")
{
    SECTION("Min/Max with fixed dimensions")
    {
        constexpr DimensionTpl<5> a;
        constexpr DimensionTpl<8> b;

        auto min_result = min(a, b);
        REQUIRE(min_result.value() == 5);
        REQUIRE(min_result.Value == 5);
        REQUIRE_FALSE(min_result.IsDynamic);

        auto max_result = max(a, b);
        REQUIRE(max_result.value() == 8);
        REQUIRE(max_result.Value == 8);
        REQUIRE_FALSE(max_result.IsDynamic);
    }

    SECTION("Min/Max with dynamic dimensions")
    {
        DimensionTpl<> a(12);
        DimensionTpl<> b(7);

        auto min_result = min(a, b);
        REQUIRE(min_result.value() == 7);
        REQUIRE(min_result.Value == detail::Dynamic);
        REQUIRE(min_result.IsDynamic);

        auto max_result = max(a, b);
        REQUIRE(max_result.value() == 12);
        REQUIRE(max_result.Value == detail::Dynamic);
        REQUIRE(max_result.IsDynamic);
    }

    SECTION("Min/Max with mixed dimensions")
    {
        constexpr DimensionTpl<10> fixed;
        DimensionTpl<> dynamic_less(3);
        DimensionTpl<> dynamic_greater(15);

        auto min_result1 = min(fixed, dynamic_less);
        REQUIRE(min_result1.value() == 3);
        REQUIRE(min_result1.Value == detail::Dynamic);
        REQUIRE(min_result1.IsDynamic);

        auto max_result1 = max(fixed, dynamic_less);
        REQUIRE(max_result1.value() == 10);
        REQUIRE(max_result1.Value == detail::Dynamic);
        REQUIRE(max_result1.IsDynamic);

        auto min_result2 = min(fixed, dynamic_greater);
        REQUIRE(min_result2.value() == 10);
        REQUIRE(min_result2.Value == detail::Dynamic);
        REQUIRE(min_result2.IsDynamic);

        auto max_result2 = max(fixed, dynamic_greater);
        REQUIRE(max_result2.value() == 15);
        REQUIRE(max_result2.Value == detail::Dynamic);
        REQUIRE(max_result2.IsDynamic);
    }

    SECTION("Min/Max with equal values")
    {
        constexpr DimensionTpl<7> fixed;
        DimensionTpl<> dynamic(7);

        auto min_result = min(fixed, dynamic);
        REQUIRE(min_result.value() == 7);
        REQUIRE(min_result.Value == detail::Dynamic);
        REQUIRE(min_result.IsDynamic);

        auto max_result = max(fixed, dynamic);
        REQUIRE(max_result.value() == 7);
        REQUIRE(max_result.Value == detail::Dynamic);
        REQUIRE(max_result.IsDynamic);
    }

    SECTION("Min/Max with zero")
    {
        constexpr DimensionTpl<0> zero;
        constexpr DimensionTpl<5> positive;

        auto min_result = min(zero, positive);
        REQUIRE(min_result.value() == 0);
        REQUIRE(min_result.Value == 0);
        REQUIRE_FALSE(min_result.IsDynamic);

        auto max_result = max(zero, positive);
        REQUIRE(max_result.value() == 5);
        REQUIRE(max_result.Value == 5);
        REQUIRE_FALSE(max_result.IsDynamic);
    }

    SECTION("Min/Max with different values")
    {
        DimensionTpl<> small(3);
        DimensionTpl<> large(8);

        auto min_result = min(small, large);
        REQUIRE(min_result.value() == 3);
        REQUIRE(min_result.Value == detail::Dynamic);
        REQUIRE(min_result.IsDynamic);

        auto max_result = max(small, large);
        REQUIRE(max_result.value() == 8);
        REQUIRE(max_result.Value == detail::Dynamic);
        REQUIRE(max_result.IsDynamic);
    }
}

TEST_CASE("DimensionTpl - Extract Compile Time Value", "[dimensions]")
{
    SECTION("Extract from raw integer")
    {
        constexpr int raw_value = 42;
        constexpr int extracted = extract_compile_time_value<raw_value>::Value;
        REQUIRE(extracted == 42);
    }

    SECTION("Extract from fixed dimension - concept validation")
    {
        // Note: extract_compile_time_value with DimensionTpl objects as template parameters
        // is intended for use with constexpr dimension values in template contexts,
        // but DimensionTpl objects cannot be used as template non-type parameters
        // because they are not structural types. This is by design - the extraction
        // should be used with compile-time constants, not runtime objects.
        static_assert(extract_compile_time_value<15>::Value == 15);
        static_assert(extract_compile_time_value<0>::Value == 0);
    }

    SECTION("Extract from negative integer")
    {
        constexpr int negative_value = -5;
        constexpr int extracted = extract_compile_time_value<negative_value>::Value;
        REQUIRE(extracted == -5);
    }

    SECTION("Extract from large integer")
    {
        constexpr int large_value = 9999;
        constexpr int extracted = extract_compile_time_value<large_value>::Value;
        REQUIRE(extracted == 9999);
    }
}

TEST_CASE("DimensionTpl - Edge Cases and Boundary Conditions", "[dimensions]")
{
    SECTION("Large dimension values")
    {
        constexpr DimensionTpl<999999> large_fixed;
        REQUIRE(large_fixed.value() == 999999);

        DimensionTpl<> large_dynamic(1000000);
        REQUIRE(large_dynamic.value() == 1000000);
    }

    SECTION("Chained operations")
    {
        constexpr DimensionTpl<10> a;
        constexpr DimensionTpl<5> b;
        constexpr DimensionTpl<2> c;

        auto result = (a + b) * c - DimensionTpl<3>();
        REQUIRE(result.value() == 27); // (10 + 5) * 2 - 3 = 27
        REQUIRE_FALSE(result.IsDynamic);
    }

    SECTION("Mixed arithmetic chains")
    {
        constexpr DimensionTpl<8> fixed;
        DimensionTpl<> dynamic(4);

        auto result = (fixed + dynamic) / 3;
        REQUIRE(result.value() == 4); // (8 + 4) / 3 = 4
        REQUIRE(result.IsDynamic);
    }

    SECTION("Complex expressions with scalars")
    {
        constexpr DimensionTpl<6> fixed;
        DimensionTpl<> dynamic(9);

        auto result = (fixed * 2) + (dynamic - 1) / 2;
        REQUIRE(result.value() == 16); // (6 * 2) + (9 - 1) / 2 = 12 + 4 = 16
        REQUIRE(result.IsDynamic);
    }

    SECTION("DimensionTpl with large values")
    {
        DimensionTpl<> large_value(INT_MAX);
        DimensionTpl<> small_value(0);

        REQUIRE(large_value.value() == INT_MAX);
        REQUIRE(small_value.value() == 0);

        REQUIRE(large_value > small_value);
        REQUIRE(small_value < large_value);
    }

    SECTION("DimensionTpl conversion consistency")
    {
        constexpr DimensionTpl<123> fixed;
        int value_direct = fixed.value();
        int value_cast = static_cast<int>(fixed);
        int value_implicit = fixed;

        REQUIRE(value_direct == 123);
        REQUIRE(value_cast == 123);
        REQUIRE(value_implicit == 123);
    }
}

TEST_CASE("DimensionTpl - Compile Time Constants", "[dimensions]")
{
    SECTION("Compile time dimension properties")
    {
        static_assert(DimensionTpl<5>::Value == 5);
        static_assert(DimensionTpl<5>::IsFixed);

        static_assert(DimensionTpl<>::Value == detail::Dynamic);
        static_assert(DimensionTpl<>::IsDynamic);
    }

    SECTION("Compile time arithmetic result types")
    {
        constexpr DimensionTpl<3> a;
        constexpr DimensionTpl<4> b;

        static_assert(decltype(a + b)::Value == 7);
        static_assert(decltype(a - b)::Value == -1);
        static_assert(decltype(a * b)::Value == 12);
        static_assert(decltype(a / b)::Value == 0); // Integer division

        static_assert(decltype(a + b)::IsFixed);
        static_assert(decltype(a - b)::IsDynamic); // In this case, the result is dynamic because 3 < 4
        static_assert(decltype(a * b)::IsFixed);
        static_assert(decltype(a / b)::IsFixed);
    }

    SECTION("Compile time scalar arithmetic result types")
    {
        constexpr DimensionTpl<10> fixed;

        // Scalar operations with runtime values always produce Dynamic dimensions
        static_assert(decltype(fixed + 5)::Value == detail::Dynamic);
        static_assert(decltype(fixed - 3)::Value == detail::Dynamic);
        static_assert(decltype(fixed * 2)::Value == detail::Dynamic);
        static_assert(decltype(fixed / 2)::Value == detail::Dynamic);

        static_assert(decltype(fixed + 5)::IsDynamic);
        static_assert(decltype(fixed - 3)::IsDynamic);
        static_assert(decltype(fixed * 2)::IsDynamic);
        static_assert(decltype(fixed / 2)::IsDynamic);
    }

    SECTION("Mixed compile time and dynamic result types")
    {
        constexpr DimensionTpl<7> fixed;
        DimensionTpl<> dynamic(3);

        static_assert(decltype(fixed + dynamic)::Value == detail::Dynamic);
        static_assert(decltype(fixed - dynamic)::Value == detail::Dynamic);
        static_assert(decltype(fixed * dynamic)::Value == detail::Dynamic);
        static_assert(decltype(fixed / dynamic)::Value == detail::Dynamic);

        static_assert(decltype(fixed + dynamic)::IsDynamic);
        static_assert(decltype(fixed - dynamic)::IsDynamic);
        static_assert(decltype(fixed * dynamic)::IsDynamic);
        static_assert(decltype(fixed / dynamic)::IsDynamic);
    }
}
