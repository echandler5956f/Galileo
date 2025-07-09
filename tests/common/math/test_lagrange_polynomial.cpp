#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/common/math/lagrange-polynomial.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <sstream>
#include <chrono>

using namespace galileo;
using namespace Catch::Matchers;

template<typename NumScalar>
constexpr NumScalar TOLERANCE = std::numeric_limits<NumScalar>::epsilon() * 100;

// Test fixture for common polynomial setups
template<typename NumScalar>
class LagrangePolynomialTestFixture
{
public:
    using PolynomialType = LagrangePolynomialTpl<NumScalar>;
    using VectorType = Eigen::Matrix<NumScalar, Eigen::Dynamic, 1>;

    // Common test polynomials
    PolynomialType constant_poly;      // $P(x) = 5$
    PolynomialType linear_poly;        // $P(x) = 2x + 3$
    PolynomialType quadratic_poly;     // $P(x) = x^2 - 4x + 3$
    PolynomialType cubic_poly;         // $P(x) = 2x^3 - 3x^2 + x - 1$
    PolynomialType zero_poly;          // $P(x) = 0$

    LagrangePolynomialTestFixture()
        : constant_poly(make_coeffs({5}))
        , linear_poly(make_coeffs({2, 3}))
        , quadratic_poly(make_coeffs({1, -4, 3}))
        , cubic_poly(make_coeffs({2, -3, 1, -1}))
        , zero_poly(make_coeffs({0}))
    {
    }

private:
    VectorType make_coeffs(std::initializer_list<NumScalar> coeffs)
    {
        VectorType vec(coeffs.size());
        std::copy(coeffs.begin(), coeffs.end(), vec.data());
        return vec;
    }
};

TEST_CASE("LagrangePolynomial - Basic Construction and Properties", "[lagrange_polynomial]")
{
    SECTION("Construction from coefficient vector")
    {
        Eigen::Vector3d coeffs;
        coeffs << 1.0, 2.0, 3.0;  // $P(x) = x^2 + 2x + 3$

        LagrangePolynomialTpl<double> poly(coeffs);
        REQUIRE(poly.get_N() == 3);
        REQUIRE(poly.get_coeffs().size() == 3);

        const auto& stored_coeffs = poly.get_coeffs();
        for (int i = 0; i < 3; ++i)
        {
            REQUIRE(stored_coeffs(i) == coeffs(i));
        }
    }

    SECTION("Construction from different vector types")
    {
        // Fixed-size vector
        Eigen::Vector4f fixed_coeffs;
        fixed_coeffs << 1.0f, -2.0f, 3.0f, -4.0f;
        LagrangePolynomialTpl<float> fixed_poly(fixed_coeffs);
        REQUIRE(fixed_poly.get_N() == 4);

        // Dynamic vector
        Eigen::VectorXd dynamic_coeffs(5);
        dynamic_coeffs << 1.0, 2.0, 3.0, 4.0, 5.0;
        LagrangePolynomialTpl<double> dynamic_poly(dynamic_coeffs);
        REQUIRE(dynamic_poly.get_N() == 5);

        // Row vector
        Eigen::RowVector3d row_coeffs;
        row_coeffs << 7.0, 8.0, 9.0;
        LagrangePolynomialTpl<double> row_poly(row_coeffs);
        REQUIRE(row_poly.get_N() == 3);
    }

    SECTION("Template parameter validation")
    {
        static_assert(std::is_same_v<LagrangePolynomialTpl<float>::NumScalar, float>);
        static_assert(LagrangePolynomialTpl<double>::Options == 0);
        static_assert(std::is_same_v<LagrangePolynomialTpl<double>::Polynomial,
                                     LagrangePolynomialTpl<double>>);
    }

    SECTION("Empty polynomial handling")
    {
        Eigen::VectorXd empty_coeffs(0);
        LagrangePolynomialTpl<double> empty_poly(empty_coeffs);

        REQUIRE(empty_poly.get_N() == 0);
        REQUIRE(empty_poly.evaluate(1.0) == 0.0);
        REQUIRE(empty_poly(5.0) == 0.0);
    }
}

TEST_CASE_METHOD((LagrangePolynomialTestFixture<double>), "LagrangePolynomial - Evaluation", "[lagrange_polynomial]")
{
    SECTION("Constant polynomial evaluation")
    {
        REQUIRE_THAT(constant_poly.evaluate(0.0), WithinAbs(5.0, TOLERANCE<double>));
        REQUIRE_THAT(constant_poly.evaluate(1.0), WithinAbs(5.0, TOLERANCE<double>));
        REQUIRE_THAT(constant_poly.evaluate(-2.5), WithinAbs(5.0, TOLERANCE<double>));
        REQUIRE_THAT(constant_poly.evaluate(100.0), WithinAbs(5.0, TOLERANCE<double>));
    }

    SECTION("Linear polynomial evaluation")
    {
        // $P(x) = 2x + 3$
        REQUIRE_THAT(linear_poly.evaluate(0.0), WithinAbs(3.0, TOLERANCE<double>));
        REQUIRE_THAT(linear_poly.evaluate(1.0), WithinAbs(5.0, TOLERANCE<double>));
        REQUIRE_THAT(linear_poly.evaluate(-1.0), WithinAbs(1.0, TOLERANCE<double>));
        REQUIRE_THAT(linear_poly.evaluate(2.5), WithinAbs(8.0, TOLERANCE<double>));
    }

    SECTION("Quadratic polynomial evaluation")
    {
        // $P(x) = x^2 - 4x + 3$
        REQUIRE_THAT(quadratic_poly.evaluate(0.0), WithinAbs(3.0, TOLERANCE<double>));
        REQUIRE_THAT(quadratic_poly.evaluate(1.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(quadratic_poly.evaluate(2.0), WithinAbs(-1.0, TOLERANCE<double>));
        REQUIRE_THAT(quadratic_poly.evaluate(3.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(quadratic_poly.evaluate(-1.0), WithinAbs(8.0, TOLERANCE<double>));
    }

    SECTION("Cubic polynomial evaluation")
    {
        // $P(x) = 2x^3 - 3x^2 + x - 1$
        REQUIRE_THAT(cubic_poly.evaluate(0.0), WithinAbs(-1.0, TOLERANCE<double>));
        REQUIRE_THAT(cubic_poly.evaluate(1.0), WithinAbs(-1.0, TOLERANCE<double>));
        REQUIRE_THAT(cubic_poly.evaluate(2.0), WithinAbs(5.0, TOLERANCE<double>));
        REQUIRE_THAT(cubic_poly.evaluate(-1.0), WithinAbs(-7.0, TOLERANCE<double>));
    }

    SECTION("Operator() evaluation")
    {
        // Should be equivalent to evaluate()
        REQUIRE_THAT(linear_poly(1.5), WithinAbs(linear_poly.evaluate(1.5), TOLERANCE<double>));
        REQUIRE_THAT(quadratic_poly(-2.0), WithinAbs(quadratic_poly.evaluate(-2.0), TOLERANCE<double>));
        REQUIRE_THAT(cubic_poly(0.5), WithinAbs(cubic_poly.evaluate(0.5), TOLERANCE<double>));
    }

    SECTION("Zero polynomial evaluation")
    {
        REQUIRE_THAT(zero_poly.evaluate(0.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(zero_poly.evaluate(1.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(zero_poly.evaluate(-100.0), WithinAbs(0.0, TOLERANCE<double>));
    }
}

TEST_CASE_METHOD((LagrangePolynomialTestFixture<double>), "LagrangePolynomial - Derivatives", "[lagrange_polynomial]")
{
    SECTION("Constant polynomial derivative")
    {
        auto derivative = constant_poly.derivative();

        REQUIRE(derivative.get_N() == 1);
        REQUIRE_THAT(derivative.evaluate(0.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(5.0), WithinAbs(0.0, TOLERANCE<double>));
    }

    SECTION("Linear polynomial derivative")
    {
        // $P(x) = 2x + 3$, $P'(x) = 2$
        auto derivative = linear_poly.derivative();

        REQUIRE(derivative.get_N() == 1);
        REQUIRE_THAT(derivative.evaluate(0.0), WithinAbs(2.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(10.0), WithinAbs(2.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(-5.0), WithinAbs(2.0, TOLERANCE<double>));
    }

    SECTION("Quadratic polynomial derivative")
    {
        // $P(x) = x^2 - 4x + 3$, $P'(x) = 2x - 4$
        auto derivative = quadratic_poly.derivative();

        REQUIRE(derivative.get_N() == 2);
        REQUIRE_THAT(derivative.evaluate(0.0), WithinAbs(-4.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(1.0), WithinAbs(-2.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(2.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(3.0), WithinAbs(2.0, TOLERANCE<double>));
    }

    SECTION("Cubic polynomial derivative")
    {
        // $P(x) = 2x^3 - 3x^2 + x - 1$, $P'(x) = 6x^2 - 6x + 1$
        auto derivative = cubic_poly.derivative();

        REQUIRE(derivative.get_N() == 3);
        REQUIRE_THAT(derivative.evaluate(0.0), WithinAbs(1.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(1.0), WithinAbs(1.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(-1.0), WithinAbs(13.0, TOLERANCE<double>));
    }

    SECTION("Second derivative")
    {
        // $P(x) = x^2 - 4x + 3$, $P'(x) = 2x - 4$, $P''(x) = 2$
        auto second_derivative = quadratic_poly.derivative().derivative();

        REQUIRE(second_derivative.get_N() == 1);
        REQUIRE_THAT(second_derivative.evaluate(0.0), WithinAbs(2.0, TOLERANCE<double>));
        REQUIRE_THAT(second_derivative.evaluate(100.0), WithinAbs(2.0, TOLERANCE<double>));
    }

    SECTION("Zero polynomial derivative")
    {
        auto derivative = zero_poly.derivative();

        REQUIRE(derivative.get_N() == 1);
        REQUIRE_THAT(derivative.evaluate(0.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(derivative.evaluate(1.0), WithinAbs(0.0, TOLERANCE<double>));
    }
}

TEST_CASE_METHOD((LagrangePolynomialTestFixture<double>), "LagrangePolynomial - Integration", "[lagrange_polynomial]")
{
    SECTION("Constant polynomial integral")
    {
        // $P(x) = 5$, $\int P(x)dx = 5x + C$
        auto integral = constant_poly.integral();

        REQUIRE(integral.get_N() == 2);
        // The integral coefficients should be [5, 0] for $5x + 0$
        REQUIRE_THAT(integral.evaluate(1.0) - integral.evaluate(0.0), WithinAbs(5.0, TOLERANCE<double>));
        REQUIRE_THAT(integral.evaluate(2.0) - integral.evaluate(0.0), WithinAbs(10.0, TOLERANCE<double>));
    }

    SECTION("Linear polynomial integral")
    {
        // $P(x) = 2x + 3$, $\int P(x)dx = x^2 + 3x + C$
        auto integral = linear_poly.integral();

        REQUIRE(integral.get_N() == 3);
        // Check by evaluating derivative
        auto derivative_of_integral = integral.derivative();
        REQUIRE_THAT(derivative_of_integral.evaluate(0.0), WithinAbs(linear_poly.evaluate(0.0), TOLERANCE<double>));
        REQUIRE_THAT(derivative_of_integral.evaluate(1.0), WithinAbs(linear_poly.evaluate(1.0), TOLERANCE<double>));
        REQUIRE_THAT(derivative_of_integral.evaluate(2.0), WithinAbs(linear_poly.evaluate(2.0), TOLERANCE<double>));
    }

    SECTION("Quadratic polynomial integral")
    {
        // $P(x) = x^2 - 4x + 3$, $\int P(x)dx = x^3/3 - 2x^2 + 3x + C$
        auto integral = quadratic_poly.integral();

        REQUIRE(integral.get_N() == 4);
        // Verify by taking derivative
        auto derivative_of_integral = integral.derivative();
        for (double x : {-2.0, -1.0, 0.0, 1.0, 2.0, 3.0})
        {
            REQUIRE_THAT(derivative_of_integral.evaluate(x),
                        WithinAbs(quadratic_poly.evaluate(x), TOLERANCE<double>));
        }
    }

    SECTION("Definite integration")
    {
        // Test $\int_0^1 (2x + 3) dx = [x^2 + 3x]_0^1 = 1 + 3 = 4$
        double result = linear_poly.integrate(0.0, 1.0);
        REQUIRE_THAT(result, WithinAbs(4.0, TOLERANCE<double>));

        // Test $\int_1^2 (x^2 - 4x + 3) dx$
        double expected = (8.0/3.0 - 8.0 + 6.0) - (1.0/3.0 - 2.0 + 3.0);
        double actual = quadratic_poly.integrate(1.0, 2.0);
        REQUIRE_THAT(actual, WithinAbs(expected, TOLERANCE<double>));
    }

    SECTION("Integration bounds")
    {
        // $\int_a^b f(x) dx = -\int_b^a f(x) dx$
        double forward = quadratic_poly.integrate(0.0, 2.0);
        double backward = quadratic_poly.integrate(2.0, 0.0);
        REQUIRE_THAT(forward, WithinAbs(-backward, TOLERANCE<double>));

        // $\int_a^a f(x) dx = 0$
        double same_bounds = cubic_poly.integrate(1.5, 1.5);
        REQUIRE_THAT(same_bounds, WithinAbs(0.0, TOLERANCE<double>));
    }
}

TEST_CASE_METHOD((LagrangePolynomialTestFixture<double>), "LagrangePolynomial - Arithmetic Operations", "[lagrange_polynomial]")
{
    SECTION("Addition of polynomials")
    {
        // $(2x + 3) + (x^2 - 4x + 3) = x^2 - 2x + 6$
        auto sum = linear_poly + quadratic_poly;

        REQUIRE(sum.get_N() == 3);
        REQUIRE_THAT(sum.evaluate(0.0), WithinAbs(6.0, TOLERANCE<double>));
        REQUIRE_THAT(sum.evaluate(1.0), WithinAbs(5.0, TOLERANCE<double>));
        REQUIRE_THAT(sum.evaluate(2.0), WithinAbs(6.0, TOLERANCE<double>));
    }

    SECTION("Addition assignment")
    {
        auto poly_copy = linear_poly;
        poly_copy += quadratic_poly;

        auto expected = linear_poly + quadratic_poly;
        for (double x : {-1.0, 0.0, 1.0, 2.0})
        {
            REQUIRE_THAT(poly_copy.evaluate(x), WithinAbs(expected.evaluate(x), TOLERANCE<double>));
        }
    }

    SECTION("Subtraction of polynomials")
    {
        // $(x^2 - 4x + 3) - (2x + 3) = x^2 - 6x$
        auto diff = quadratic_poly - linear_poly;

        REQUIRE(diff.get_N() == 3);
        REQUIRE_THAT(diff.evaluate(0.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(diff.evaluate(1.0), WithinAbs(-5.0, TOLERANCE<double>));
        REQUIRE_THAT(diff.evaluate(2.0), WithinAbs(-8.0, TOLERANCE<double>));
    }

    SECTION("Subtraction assignment")
    {
        auto poly_copy = quadratic_poly;
        poly_copy -= linear_poly;

        auto expected = quadratic_poly - linear_poly;
        for (double x : {-1.0, 0.0, 1.0, 2.0})
        {
            REQUIRE_THAT(poly_copy.evaluate(x), WithinAbs(expected.evaluate(x), TOLERANCE<double>));
        }
    }

        SECTION("Multiplication of polynomials")
    {
        // $(2x + 3) \times (x^2 - 4x + 3) = 2x^3 + 3x^2 - 8x^2 - 12x + 6x + 9$
        //                           $= 2x^3 - 5x^2 - 6x + 9$
        auto product = linear_poly * quadratic_poly;

        REQUIRE(product.get_N() == 4);
        REQUIRE_THAT(product.evaluate(0.0), WithinAbs(9.0, TOLERANCE<double>));
        REQUIRE_THAT(product.evaluate(1.0), WithinAbs(0.0, TOLERANCE<double>));
        // At $x = -1$: $2(-1)^3 - 5(-1)^2 - 6(-1) + 9 = -2 - 5 + 6 + 9 = 8$
        REQUIRE_THAT(product.evaluate(-1.0), WithinAbs(8.0, TOLERANCE<double>));
    }

    SECTION("Multiplication assignment")
    {
        auto poly_copy = linear_poly;
        poly_copy *= quadratic_poly;

        auto expected = linear_poly * quadratic_poly;
        for (double x : {-1.0, 0.0, 1.0, 2.0})
        {
            REQUIRE_THAT(poly_copy.evaluate(x), WithinAbs(expected.evaluate(x), TOLERANCE<double>));
        }
    }

    SECTION("Addition with different sizes")
    {
        // constant + cubic should work
        auto result = constant_poly + cubic_poly;
        REQUIRE(result.get_N() == 4);

        for (double x : {-1.0, 0.0, 1.0, 2.0})
        {
            double expected = constant_poly.evaluate(x) + cubic_poly.evaluate(x);
            REQUIRE_THAT(result.evaluate(x), WithinAbs(expected, TOLERANCE<double>));
        }
    }

    SECTION("Operations with zero polynomial")
    {
        auto sum_with_zero = linear_poly + zero_poly;
        auto zero_sum = zero_poly + linear_poly;

        for (double x : {-2.0, 0.0, 1.0, 5.0})
        {
            REQUIRE_THAT(sum_with_zero.evaluate(x), WithinAbs(linear_poly.evaluate(x), TOLERANCE<double>));
            REQUIRE_THAT(zero_sum.evaluate(x), WithinAbs(linear_poly.evaluate(x), TOLERANCE<double>));
        }
    }
}

TEST_CASE("LagrangePolynomial - Edge Cases and Special Values", "[lagrange_polynomial][edge_cases]")
{
    SECTION("Single coefficient polynomial")
    {
        Eigen::VectorXd single_coeff(1);
        single_coeff << 42.0;
        LagrangePolynomialTpl<double> single_poly(single_coeff);

        REQUIRE(single_poly.get_N() == 1);
        REQUIRE_THAT(single_poly.evaluate(0.0), WithinAbs(42.0, TOLERANCE<double>));
        REQUIRE_THAT(single_poly.evaluate(100.0), WithinAbs(42.0, TOLERANCE<double>));

        auto derivative = single_poly.derivative();
        REQUIRE_THAT(derivative.evaluate(0.0), WithinAbs(0.0, TOLERANCE<double>));
    }

    SECTION("Very large coefficients")
    {
        Eigen::Vector3d large_coeffs;
        large_coeffs << 1e10, -1e10, 1e10;
        LagrangePolynomialTpl<double> large_poly(large_coeffs);

        // Should handle large numbers without overflow
        REQUIRE(std::isfinite(large_poly.evaluate(0.1)));
        REQUIRE(std::isfinite(large_poly.evaluate(0.0)));
    }

    SECTION("Very small coefficients")
    {
        Eigen::Vector3d small_coeffs;
        small_coeffs << 1e-15, 1e-14, 1e-13;
        LagrangePolynomialTpl<double> small_poly(small_coeffs);

        // Should handle small numbers without underflow
        REQUIRE_THAT(small_poly.evaluate(1.0), WithinAbs(1.111e-13, TOLERANCE<double>));
    }

    SECTION("Evaluation at extreme values")
    {
        Eigen::Vector3d coeffs;
        coeffs << 1.0, 0.0, 1.0;  // $P(x) = x^2 + 1$
        LagrangePolynomialTpl<double> poly(coeffs);

        // Large positive values
        REQUIRE(poly.evaluate(1000.0) > 0.0);
        REQUIRE(std::isfinite(poly.evaluate(1000.0)));

        // Large negative values
        REQUIRE(poly.evaluate(-1000.0) > 0.0);
        REQUIRE(std::isfinite(poly.evaluate(-1000.0)));
    }
}

TEST_CASE("LagrangePolynomial - Numerical Precision and Stability", "[lagrange_polynomial][numerical]")
{
        SECTION("Horner's method precision")
    {
        // Test that Horner's method gives accurate results
        Eigen::VectorXd coeffs(10);
        for (int i = 0; i < 10; ++i)
        {
            coeffs(i) = std::pow(-1.0, i) / (i + 1.0);
        }

        LagrangePolynomialTpl<double> poly(coeffs);

        // Evaluate at a point where we can compute the expected result
        double x = 0.5;
        double expected = 0.0;

        // Manual computation: coeffs[i] corresponds to coefficient of x^(N-1-i)
        // For our 10-coefficient polynomial: coeffs[0]*x^9 + coeffs[1]*x^8 + ... + coeffs[9]*x^0
        for (int i = 0; i < 10; ++i)
        {
            int power = 9 - i;  // Power of x for coeffs[i]
            expected += coeffs(i) * std::pow(x, power);
        }

        REQUIRE_THAT(poly.evaluate(x), WithinAbs(expected, TOLERANCE<double>));
    }

    SECTION("Consistency of operations")
    {
        Eigen::Vector4d coeffs1, coeffs2;
        coeffs1 << 1.0, 2.0, 3.0, 4.0;
        coeffs2 << 5.0, 6.0, 7.0, 8.0;

        LagrangePolynomialTpl<double> p1(coeffs1);
        LagrangePolynomialTpl<double> p2(coeffs2);

        // (p1 + p2) - p1 should equal p2
        auto sum = p1 + p2;
        auto diff = sum - p1;

        for (double x : {-2.0, -0.5, 0.0, 1.0, 2.5})
        {
            REQUIRE_THAT(diff.evaluate(x), WithinAbs(p2.evaluate(x), TOLERANCE<double>));
        }
    }

    SECTION("Derivative-integral consistency")
    {
        Eigen::Vector4d coeffs;
        coeffs << 2.0, -3.0, 1.0, 5.0;
        LagrangePolynomialTpl<double> poly(coeffs);

        // $\frac{d}{dx} \int f(x)dx$ should equal $f(x)$
        auto integral = poly.integral();
        auto derivative_of_integral = integral.derivative();

        for (double x : {-1.0, 0.0, 1.0, 2.0})
        {
            REQUIRE_THAT(derivative_of_integral.evaluate(x), WithinAbs(poly.evaluate(x), TOLERANCE<double>));
        }
    }
}

TEMPLATE_TEST_CASE("LagrangePolynomial - Different Scalar Types", "[lagrange_polynomial][scalar_types]",
                   float, double, long double)
{
    using Scalar = TestType;
    constexpr Scalar tolerance = TOLERANCE<Scalar>;

    // Create a test polynomial: $P(x) = 2x^2 - 3x + 1$
    Eigen::Matrix<Scalar, 3, 1> coeffs;
    coeffs << Scalar(2), Scalar(-3), Scalar(1);

    LagrangePolynomialTpl<Scalar> poly(coeffs);

    SECTION("Basic evaluation")
    {
        REQUIRE(poly.get_N() == 3);

        // $P(0) = 1$
        REQUIRE_THAT(static_cast<double>(poly.evaluate(Scalar(0))), WithinAbs(1.0, tolerance));

        // $P(1) = 0$
        REQUIRE_THAT(static_cast<double>(poly.evaluate(Scalar(1))), WithinAbs(0.0, tolerance));

        // $P(2) = 3$
        REQUIRE_THAT(static_cast<double>(poly.evaluate(Scalar(2))), WithinAbs(3.0, tolerance));
    }

    SECTION("Derivative")
    {
        // $P'(x) = 4x - 3$
        auto derivative = poly.derivative();

        REQUIRE_THAT(static_cast<double>(derivative.evaluate(Scalar(0))), WithinAbs(-3.0, tolerance));
        REQUIRE_THAT(static_cast<double>(derivative.evaluate(Scalar(1))), WithinAbs(1.0, tolerance));
    }

    SECTION("Integration")
    {
        // $\int_0^1 (2x^2 - 3x + 1) dx = [2x^3/3 - 3x^2/2 + x]_0^1 = 2/3 - 3/2 + 1 = 1/6$
        Scalar result = poly.integrate(Scalar(0), Scalar(1));
        REQUIRE_THAT(static_cast<double>(result), WithinAbs(1.0/6.0, tolerance));
    }
}

TEST_CASE("LagrangePolynomial - Copy and Assignment", "[lagrange_polynomial][copy]")
{
    SECTION("Copy construction")
    {
        Eigen::Vector3d coeffs;
        coeffs << 1.0, -2.0, 3.0;
        LagrangePolynomialTpl<double> original(coeffs);

        LagrangePolynomialTpl<double> copy(original);

        REQUIRE(copy.get_N() == original.get_N());

        for (double x : {-1.0, 0.0, 1.0, 2.0})
        {
            REQUIRE_THAT(copy.evaluate(x), WithinAbs(original.evaluate(x), TOLERANCE<double>));
        }
    }

    SECTION("Assignment operator")
    {
        Eigen::Vector2d coeffs1, coeffs2;
        coeffs1 << 1.0, 2.0;
        coeffs2 << 3.0, 4.0;

        LagrangePolynomialTpl<double> poly1(coeffs1);
        LagrangePolynomialTpl<double> poly2(coeffs2);

        poly1 = poly2;

        REQUIRE(poly1.get_N() == poly2.get_N());

        for (double x : {-1.0, 0.0, 1.0, 2.0})
        {
            REQUIRE_THAT(poly1.evaluate(x), WithinAbs(poly2.evaluate(x), TOLERANCE<double>));
        }
    }
}

TEST_CASE("LagrangePolynomial - Stream Output", "[lagrange_polynomial][output]")
{
    SECTION("Stream output formatting")
    {
        Eigen::Vector3d coeffs;
        coeffs << 1.0, -2.0, 3.0;
        LagrangePolynomialTpl<double> poly(coeffs);

        std::ostringstream oss;
        oss << poly;

        std::string output = oss.str();

        // Should contain basic information about the polynomial
        REQUIRE(output.find("LagrangePolynomial") != std::string::npos);
        REQUIRE(output.find("degree=2") != std::string::npos);
        REQUIRE(output.find("coefficients") != std::string::npos);
    }
}

TEST_CASE("LagrangePolynomial - Performance and Large Polynomials", "[lagrange_polynomial][performance]")
{
    SECTION("Large polynomial evaluation")
    {
        constexpr int LARGE_N = 1000;
        Eigen::VectorXd large_coeffs(LARGE_N);

        // Create a polynomial with many terms
        for (int i = 0; i < LARGE_N; ++i)
        {
            large_coeffs(i) = 1.0 / (i + 1.0);
        }

        LagrangePolynomialTpl<double> large_poly(large_coeffs);

        REQUIRE(large_poly.get_N() == LARGE_N);

        // Should be able to evaluate efficiently
        auto start = std::chrono::high_resolution_clock::now();
        double result = large_poly.evaluate(0.5);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        // Should complete very quickly (less than 1ms for N=1000)
        REQUIRE(duration.count() < 1000);
        REQUIRE(std::isfinite(result));
    }

    SECTION("Large polynomial operations")
    {
        constexpr int N = 100;
        Eigen::VectorXd coeffs1(N), coeffs2(N);

        for (int i = 0; i < N; ++i)
        {
            coeffs1(i) = std::sin(i * 0.1);
            coeffs2(i) = std::cos(i * 0.1);
        }

        LagrangePolynomialTpl<double> poly1(coeffs1);
        LagrangePolynomialTpl<double> poly2(coeffs2);

        // Operations should complete efficiently
        auto sum = poly1 + poly2;
        auto product = poly1 * poly2;
        auto derivative = poly1.derivative();
        auto integral = poly1.integral();

        REQUIRE(sum.get_N() == N);
        REQUIRE(product.get_N() == 2 * N - 1);
        REQUIRE(derivative.get_N() == N - 1);
        REQUIRE(integral.get_N() == N + 1);
    }
}

TEST_CASE("LagrangePolynomial - Mathematical Properties", "[lagrange_polynomial][mathematical]")
{
    SECTION("Fundamental theorem of calculus")
    {
        // For any polynomial $P$, $\int_a^b P'(x) dx = P(b) - P(a)$
        Eigen::Vector4d coeffs;
        coeffs << 1.0, -2.0, 3.0, -1.0;  // $P(x) = x^3 - 2x^2 + 3x - 1$
        LagrangePolynomialTpl<double> poly(coeffs);

        double a = -1.0, b = 2.0;
        auto derivative = poly.derivative();
        double integral_of_derivative = derivative.integrate(a, b);
        double expected = poly.evaluate(b) - poly.evaluate(a);

        REQUIRE_THAT(integral_of_derivative, WithinAbs(expected, TOLERANCE<double>));
    }

    SECTION("Linearity of operations")
    {
        Eigen::Vector3d coeffs1, coeffs2;
        coeffs1 << 1.0, 2.0, 3.0;
        coeffs2 << 4.0, 5.0, 6.0;

        LagrangePolynomialTpl<double> p1(coeffs1);
        LagrangePolynomialTpl<double> p2(coeffs2);

        // $(p1 + p2)' = p1' + p2'$
        auto sum_derivative = (p1 + p2).derivative();
        auto derivative_sum = p1.derivative() + p2.derivative();

        for (double x : {-1.0, 0.0, 1.0, 2.0})
        {
            REQUIRE_THAT(sum_derivative.evaluate(x), WithinAbs(derivative_sum.evaluate(x), TOLERANCE<double>));
        }
    }

    SECTION("Polynomial roots and evaluation")
    {
        // $P(x) = (x-1)(x-2) = x^2 - 3x + 2$
        Eigen::Vector3d coeffs;
        coeffs << 1.0, -3.0, 2.0;
        LagrangePolynomialTpl<double> poly(coeffs);

        // Should have roots at $x = 1$ and $x = 2$
        REQUIRE_THAT(poly.evaluate(1.0), WithinAbs(0.0, TOLERANCE<double>));
        REQUIRE_THAT(poly.evaluate(2.0), WithinAbs(0.0, TOLERANCE<double>));

        // Should be positive for $x > 2$ and $x < 1$
        REQUIRE(poly.evaluate(3.0) > 0.0);
        REQUIRE(poly.evaluate(0.0) > 0.0);

        // Should be negative for $1 < x < 2$
        REQUIRE(poly.evaluate(1.5) < 0.0);
    }
}
