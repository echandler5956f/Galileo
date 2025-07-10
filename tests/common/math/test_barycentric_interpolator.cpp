#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/common/math/barycentric-interpolator.hpp"
#include "galileo/common/math/jacobi-roots.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

using namespace galileo;
using namespace Catch::Matchers;

template <typename NumScalar>
constexpr NumScalar TOLERANCE = std::numeric_limits<NumScalar>::epsilon() * 100;

// Test fixture for common interpolator setups
template <typename NumScalar, int N>
class BarycentricInterpolatorTestFixture
{
public:
    using InterpolatorType = BarycentricInterpolatorTpl<NumScalar, N>;
    using VectorType = typename InterpolatorType::VectorN;
    using MatrixType = typename InterpolatorType::MatrixN;

    // Different node distributions
    VectorType uniform_nodes;
    VectorType chebyshev_nodes;
    VectorType legendre_nodes;
    VectorType jacobi_nodes;

    // Interpolators
    std::unique_ptr<InterpolatorType> uniform_interp;
    std::unique_ptr<InterpolatorType> chebyshev_interp;
    std::unique_ptr<InterpolatorType> legendre_interp;
    std::unique_ptr<InterpolatorType> jacobi_interp;

    BarycentricInterpolatorTestFixture()
    {
        const int n = (N == Eigen::Dynamic) ? 5 : N;

        // Uniform nodes
        uniform_nodes.resize(n);
        for (int i = 0; i < n; ++i)
        {
            uniform_nodes(i) = static_cast<NumScalar>(i) / static_cast<NumScalar>(n - 1);
        }

        // Chebyshev nodes (mapped to [0, 1])
        chebyshev_nodes.resize(n);
        for (int i = 0; i < n; ++i)
        {
            NumScalar theta = M_PI * (2 * i + 1) / (2 * n);
            chebyshev_nodes(i) = 0.5 * (1.0 - std::cos(theta));
        }

        // Legendre nodes using Jacobi roots with α=β=0
        JacobiRootsTpl<NumScalar, N> legendre_roots(0.0, 0.0);
        legendre_roots.compute_roots();
        legendre_nodes = legendre_roots.get_roots();

        // General Jacobi nodes with α=0.5, β=1.5
        JacobiRootsTpl<NumScalar, N> jacobi_roots(0.5, 1.5);
        jacobi_roots.compute_roots();
        jacobi_nodes = jacobi_roots.get_roots();

        // Create interpolators
        uniform_interp = std::make_unique<InterpolatorType>(uniform_nodes);
        chebyshev_interp = std::make_unique<InterpolatorType>(chebyshev_nodes);
        legendre_interp = std::make_unique<InterpolatorType>(legendre_nodes);
        jacobi_interp = std::make_unique<InterpolatorType>(jacobi_nodes);
    }
};

struct LegendreParams
{
    static constexpr double ALPHA = 0.0;
    static constexpr double BETA = 0.0;
};

struct Chebyshev1Params
{
    static constexpr double ALPHA = -0.5;
    static constexpr double BETA = -0.5;
};

struct Chebyshev2Params
{
    static constexpr double ALPHA = 0.5;
    static constexpr double BETA = 0.5;
};

TEST_CASE("BarycentricInterpolator - Basic Construction and Properties", "[barycentric][construction]")
{
    SECTION("Fixed dimension construction")
    {
        constexpr int N = 4;
        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.33, 0.67, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        REQUIRE(interp.get_n() == N);
        REQUIRE(interp.NDim().value() == N);
        REQUIRE(interp.NDim().IsFixed);
        REQUIRE_FALSE(interp.NDim().IsDynamic);
    }

    SECTION("Dynamic dimension construction")
    {
        const int n = 6;
        Eigen::GMatrix<double, Eigen::Dynamic, 1> nodes(n);
        for (int i = 0; i < n; ++i)
        {
            nodes(i) = static_cast<double>(i) / (n - 1);
        }

        BarycentricInterpolatorTpl<double, Eigen::Dynamic> interp(nodes);

        REQUIRE(interp.get_n() == n);
        REQUIRE(interp.NDim().value() == n);
        REQUIRE_FALSE(interp.NDim().IsFixed);
        REQUIRE(interp.NDim().IsDynamic);
    }

    SECTION("Construction with different node counts")
    {
        for (int n : {1, 2, 3, 5, 10, 20})
        {
            Eigen::VectorXd nodes = Eigen::VectorXd::LinSpaced(n, 0.0, 1.0);
            BarycentricInterpolatorTpl<double, Eigen::Dynamic> interp(nodes);
            REQUIRE(interp.get_n() == n);
        }
    }

    SECTION("Template parameter validation")
    {
        using Interp1 = BarycentricInterpolatorTpl<float, 5>;
        using Interp2 = BarycentricInterpolatorTpl<double, 8, Eigen::RowMajor>;

        static_assert(Interp1::N == 5);
        static_assert(std::is_same_v<Interp1::NumScalar, float>);
        static_assert(Interp1::Options == 0);

        static_assert(Interp2::N == 8);
        static_assert(std::is_same_v<Interp2::NumScalar, double>);
        static_assert(Interp2::Options == Eigen::RowMajor);
    }
}

TEST_CASE_METHOD((BarycentricInterpolatorTestFixture<double, 5>),
                 "BarycentricInterpolator - Scalar Interpolation",
                 "[barycentric][scalar]")
{
    SECTION("Constant function interpolation")
    {
        const double constant = 3.14159;
        Eigen::GMatrix<double, 1, 5> values; // 1×5 matrix (1 output dimension, 5 nodes)
        values.setConstant(constant);

        // Test at various points
        for (double t : {0.0, 0.1, 0.5, 0.9, 1.0})
        {
            Eigen::GMatrix<double, 1, 1> result; // 1×1 matrix (scalar output)

            uniform_interp->calc(t, values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(constant, TOLERANCE<double>));

            chebyshev_interp->calc(t, values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(constant, TOLERANCE<double>));

            legendre_interp->calc(t, values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(constant, TOLERANCE<double>));
        }
    }

    SECTION("Linear function interpolation")
    {
        auto linear_func = [](double x)
        { return 2.0 * x + 1.0; };

        // Set up function values at nodes (1×5 matrices)
        Eigen::GMatrix<double, 1, 5> uniform_values, chebyshev_values, legendre_values;
        for (int i = 0; i < 5; ++i)
        {
            uniform_values(0, i) = linear_func(uniform_nodes(i));
            chebyshev_values(0, i) = linear_func(chebyshev_nodes(i));
            legendre_values(0, i) = linear_func(legendre_nodes(i));
        }

        // Test interpolation accuracy
        for (double t : {0.0, 0.25, 0.5, 0.75, 1.0})
        {
            double expected = linear_func(t);
            Eigen::GMatrix<double, 1, 1> result;

            uniform_interp->calc(t, uniform_values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(expected, TOLERANCE<double> * 10));

            chebyshev_interp->calc(t, chebyshev_values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(expected, TOLERANCE<double> * 10));

            legendre_interp->calc(t, legendre_values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(expected, TOLERANCE<double> * 10));
        }
    }

    SECTION("Quadratic function interpolation")
    {
        auto quad_func = [](double x)
        { return x * x - 2.0 * x + 3.0; };

        Eigen::GMatrix<double, 1, 5> values;
        for (int i = 0; i < 5; ++i)
        {
            values(0, i) = quad_func(chebyshev_nodes(i));
        }

        // Should be exact for polynomial of degree < N
        for (double t : {0.05, 0.15, 0.35, 0.65, 0.85, 0.95})
        {
            double expected = quad_func(t);
            Eigen::GMatrix<double, 1, 1> result;

            chebyshev_interp->calc(t, values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(expected, TOLERANCE<double> * 100));
        }
    }

    SECTION("High-degree polynomial interpolation")
    {
        // Test polynomial of degree N-1 (should be exact)
        auto poly_func = [](double x)
        {
            return 1.0 + 2.0 * x - 3.0 * x * x + 4.0 * x * x * x - 5.0 * x * x * x * x;
        };

        Eigen::GMatrix<double, 1, 5> values;
        for (int i = 0; i < 5; ++i)
        {
            values(0, i) = poly_func(legendre_nodes(i));
        }

        // Test at many points
        for (int i = 0; i <= 20; ++i)
        {
            double t = i / 20.0;
            double expected = poly_func(t);
            Eigen::GMatrix<double, 1, 1> result;

            legendre_interp->calc(t, values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(expected, TOLERANCE<double> * 1000));
        }
    }
}

TEST_CASE_METHOD((BarycentricInterpolatorTestFixture<double, 12>),
                 "BarycentricInterpolator - Transcendental Functions",
                 "[barycentric][transcendental]")
{
    SECTION("Sine function interpolation")
    {
        auto func = [](double x)
        { return std::sin(2.0 * M_PI * x); };

        Eigen::GMatrix<double, 1, 12> values; // 1×12 matrix
        for (int i = 0; i < 12; ++i)
        {
            values(0, i) = func(chebyshev_nodes(i));
        }

        // Test interpolation error
        double max_error = 0.0;
        for (int i = 0; i <= 100; ++i)
        {
            double t = i / 100.0;
            double expected = func(t);
            Eigen::GMatrix<double, 1, 1> result;

            chebyshev_interp->calc(t, values, result);
            double error = std::abs(result(0, 0) - expected);
            max_error = std::max(max_error, error);
        }

        // With 12 Chebyshev nodes, should achieve reasonable accuracy
        REQUIRE(max_error < 1e-3);
    }

    SECTION("Exponential function interpolation")
    {
        auto func = [](double x)
        { return std::exp(x); };

        Eigen::GMatrix<double, 1, 12> values;
        for (int i = 0; i < 12; ++i)
        {
            values(0, i) = func(legendre_nodes(i));
        }

        // Test at various points
        for (double t : {0.1, 0.3, 0.5, 0.7, 0.9})
        {
            double expected = func(t);
            Eigen::GMatrix<double, 1, 1> result;

            legendre_interp->calc(t, values, result);
            double relative_error = std::abs((result(0, 0) - expected) / expected);
            REQUIRE(relative_error < 1e-4);
        }
    }

    SECTION("Runge function interpolation")
    {
        // Runge function: notorious for oscillations with uniform nodes
        auto runge = [](double x)
        {
            double shifted = 2.0 * x - 1.0; // Map [0,1] to [-1,1]
            return 1.0 / (1.0 + 25.0 * shifted * shifted);
        };

        // Compare uniform vs Chebyshev nodes
        Eigen::GMatrix<double, 1, 12> uniform_values, chebyshev_values;
        for (int i = 0; i < 12; ++i)
        {
            uniform_values(0, i) = runge(uniform_nodes(i));
            chebyshev_values(0, i) = runge(chebyshev_nodes(i));
        }

        // Compute max errors
        double uniform_max_error = 0.0;
        double chebyshev_max_error = 0.0;

        for (int i = 0; i <= 100; ++i)
        {
            double t = i / 100.0;
            double expected = runge(t);
            Eigen::GMatrix<double, 1, 1> result;

            uniform_interp->calc(t, uniform_values, result);
            uniform_max_error = std::max(uniform_max_error, std::abs(result(0, 0) - expected));

            chebyshev_interp->calc(t, chebyshev_values, result);
            chebyshev_max_error = std::max(chebyshev_max_error, std::abs(result(0, 0) - expected));
        }

        // Chebyshev nodes should perform significantly better
        REQUIRE(chebyshev_max_error < uniform_max_error);
    }
}

TEST_CASE("BarycentricInterpolator - Edge Cases", "[barycentric][edge_cases]")
{
    SECTION("Interpolation at nodes")
    {
        constexpr int N = 4;
        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.3, 0.7, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        Eigen::GMatrix<double, 1, N> values; // 1×4 matrix
        values << 1.0, -2.0, 3.5, -0.5;

        // Interpolation at nodes should return exact values
        for (int i = 0; i < N; ++i)
        {
            Eigen::GMatrix<double, 1, 1> result;
            interp.calc(nodes(i), values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(values(0, i), TOLERANCE<double>));
        }
    }

    SECTION("Interpolation near nodes")
    {
        constexpr int N = 3;
        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.5, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        Eigen::GMatrix<double, 1, N> values;
        values << 2.0, -1.0, 3.0;

        // Test interpolation very close to nodes
        for (int i = 0; i < N; ++i)
        {
            double eps = std::numeric_limits<double>::epsilon() * 10;

            // Just before node
            if (nodes(i) - eps >= 0.0)
            {
                Eigen::GMatrix<double, 1, 1> result;
                interp.calc(nodes(i) - eps, values, result);
                REQUIRE_THAT(result(0, 0), WithinAbs(values(0, i), 1e-6));
            }

            // Just after node
            if (nodes(i) + eps <= 1.0)
            {
                Eigen::GMatrix<double, 1, 1> result;
                interp.calc(nodes(i) + eps, values, result);
                REQUIRE_THAT(result(0, 0), WithinAbs(values(0, i), 1e-6));
            }
        }
    }

    SECTION("Boundary interpolation (t=0 and t=1)")
    {
        Eigen::GMatrix<double, 5, 1> nodes = Eigen::GMatrix<double, 5, 1>::LinSpaced(5, 0.0, 1.0);
        BarycentricInterpolatorTpl<double, 5> interp(nodes);

        Eigen::GMatrix<double, 1, 5> values;
        values << 1.5, 2.0, 0.5, -1.0, 3.5;

        Eigen::GMatrix<double, 1, 1> result;

        // t = 0
        interp.calc(0.0, values, result);
        REQUIRE_THAT(result(0, 0), WithinAbs(values(0, 0), TOLERANCE<double>));

        // t = 1
        interp.calc(1.0, values, result);
        REQUIRE_THAT(result(0, 0), WithinAbs(values(0, 4), TOLERANCE<double>));
    }

    SECTION("Single node interpolation (N=1)")
    {
        Eigen::GMatrix<double, 1, 1> nodes;
        nodes << 0.5;

        BarycentricInterpolatorTpl<double, 1> interp(nodes);

        Eigen::GMatrix<double, 1, 1> values;
        values << 2.718;

        // Should return constant value for any t
        for (double t : {0.0, 0.25, 0.5, 0.75, 1.0})
        {
            Eigen::GMatrix<double, 1, 1> result;
            interp.calc(t, values, result);
            REQUIRE_THAT(result(0, 0), WithinAbs(values(0, 0), TOLERANCE<double>));
        }
    }

    SECTION("Two node interpolation (N=2)")
    {
        Eigen::GMatrix<double, 2, 1> nodes;
        nodes << 0.0, 1.0;

        BarycentricInterpolatorTpl<double, 2> interp(nodes);

        Eigen::GMatrix<double, 1, 2> values;
        values << 1.0, 3.0;

        // Should perform linear interpolation
        for (double t : {0.0, 0.25, 0.5, 0.75, 1.0})
        {
            Eigen::GMatrix<double, 1, 1> result;
            interp.calc(t, values, result);
            double expected = 1.0 + 2.0 * t;
            REQUIRE_THAT(result(0, 0), WithinAbs(expected, TOLERANCE<double>));
        }
    }
}

TEST_CASE("BarycentricInterpolator - Multi-dimensional Interpolation", "[barycentric][multidim]")
{
    SECTION("Vector-valued function interpolation")
    {
        constexpr int N = 12;
        constexpr int DIM = 3;

        Eigen::GMatrix<double, N, 1> nodes;
        for (int i = 0; i < N; ++i)
        {
            double theta = M_PI * (2 * i + 1) / (2 * N);
            nodes(i) = 0.5 * (1.0 - std::cos(theta));
        }

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Create vector-valued data (3×12 matrix)
        Eigen::GMatrix<double, DIM, N> values;
        for (int i = 0; i < N; ++i)
        {
            values(0, i) = std::sin(2 * M_PI * nodes(i));
            values(1, i) = std::cos(2 * M_PI * nodes(i));
            values(2, i) = nodes(i) * nodes(i);
        }

        // Test interpolation
        for (double t : {0.1, 0.3, 0.5, 0.7, 0.9})
        {
            Eigen::GMatrix<double, DIM, 1> result;
            interp.calc(t, values, result);

            // Check each dimension
            REQUIRE_THAT(result(0, 0), WithinAbs(std::sin(2 * M_PI * t), 1e-3));
            REQUIRE_THAT(result(1, 0), WithinAbs(std::cos(2 * M_PI * t), 1e-3));
            REQUIRE_THAT(result(2, 0), WithinAbs(t * t, 1e-6)); // Polynomial should be very accurate
        }
    }

    SECTION("Matrix-valued interpolation")
    {
        constexpr int N = 4;

        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.33, 0.67, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Create matrix data (4×4 matrix, representing flattened 2x2 matrices at each node)
        Eigen::GMatrix<double, 4, N> values;
        for (int i = 0; i < N; ++i)
        {
            double t = nodes(i);
            values(0, i) = t;       // (0,0) element
            values(1, i) = t * t;   // (0,1) element
            values(2, i) = 1.0 - t; // (1,0) element
            values(3, i) = 2.0 * t; // (1,1) element
        }

        // Test interpolation
        double test_t = 0.5;
        Eigen::GMatrix<double, 4, 1> result;
        interp.calc(test_t, values, result);

        // Verify results
        REQUIRE_THAT(result(0, 0), WithinAbs(test_t, 1e-3));
        REQUIRE_THAT(result(1, 0), WithinAbs(test_t * test_t, 1e-3));
        REQUIRE_THAT(result(2, 0), WithinAbs(1.0 - test_t, 1e-3));
        REQUIRE_THAT(result(3, 0), WithinAbs(2.0 * test_t, 1e-3));
    }
}

TEST_CASE("BarycentricInterpolator - calcDiff Method", "[barycentric][derivatives]")
{
    SECTION("Basic sensitivity computation")
    {
        constexpr int N = 3;
        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.5, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        Eigen::GMatrix<double, 1, N> values; // 1×3 matrix
        values << 1.0, 2.0, 3.0;

        double t = 0.3;
        Eigen::GMatrix<double, 1, 1 * N> du_dw; // 1×3 matrix
        interp.calcDiff(t, values, du_dw);

        // Check dimensions
        REQUIRE(du_dw.rows() == 1);
        REQUIRE(du_dw.cols() == 3);

        // Sum of sensitivities should be 1 (partition of unity)
        double sum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            sum += du_dw(0, i);
        }
        REQUIRE_THAT(sum, WithinAbs(1.0, TOLERANCE<double>));
    }

    SECTION("Sensitivity at interpolation nodes")
    {
        constexpr int N = 4;
        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.3, 0.7, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        Eigen::GMatrix<double, 2, N> values; // 2×4 matrix
        values.row(0) << 1.0, 2.0, 3.0, 4.0;
        values.row(1) << -1.0, 0.0, 1.0, 2.0;

        // Test at each node
        for (int k = 0; k < N; ++k)
        {
            Eigen::GMatrix<double, 2, 2 * N> du_dw; // 2×8 matrix
            interp.calcDiff(nodes(k), values, du_dw);

            // At node k, the sensitivity should be identity for column k
            // The output is organized as blocks of 2×2 (one for each node)
            for (int j = 0; j < N; ++j)
            {
                Eigen::Matrix2d block = du_dw.block<2, 2>(0, j * 2);

                if (j == k)
                {
                    // Should be identity matrix
                    REQUIRE_THAT(block(0, 0), WithinAbs(1.0, TOLERANCE<double>));
                    REQUIRE_THAT(block(1, 1), WithinAbs(1.0, TOLERANCE<double>));
                    REQUIRE_THAT(block(0, 1), WithinAbs(0.0, TOLERANCE<double>));
                    REQUIRE_THAT(block(1, 0), WithinAbs(0.0, TOLERANCE<double>));
                }
                else
                {
                    // Should be zero matrix
                    REQUIRE_THAT(block(0, 0), WithinAbs(0.0, TOLERANCE<double>));
                    REQUIRE_THAT(block(1, 1), WithinAbs(0.0, TOLERANCE<double>));
                    REQUIRE_THAT(block(0, 1), WithinAbs(0.0, TOLERANCE<double>));
                    REQUIRE_THAT(block(1, 0), WithinAbs(0.0, TOLERANCE<double>));
                }
            }
        }
    }

    SECTION("Finite difference verification")
    {
        constexpr int N = 3;
        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.5, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        Eigen::GMatrix<double, 2, N> values; // 2×3 matrix
        values.row(0) << 1.0, 2.5, 0.5;
        values.row(1) << -1.0, 1.0, 2.0;

        double t = 0.6;
        const double eps = 1e-7;

        // Compute analytical derivatives
        Eigen::GMatrix<double, 2, 2 * N> du_dw; // 2×6 matrix
        interp.calcDiff(t, values, du_dw);

        // Compute finite difference approximation
        Eigen::GMatrix<double, 2, 1> u0, u_plus, u_minus;
        interp.calc(t, values, u0);

        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < 2; ++i)
            {
                // Perturb values
                Eigen::GMatrix<double, 2, N> values_plus = values;
                Eigen::GMatrix<double, 2, N> values_minus = values;
                values_plus(i, j) += eps;
                values_minus(i, j) -= eps;

                interp.calc(t, values_plus, u_plus);
                interp.calc(t, values_minus, u_minus);

                double fd_derivative = (u_plus(i, 0) - u_minus(i, 0)) / (2 * eps);

                // The analytical derivative is in the diagonal of the j-th block
                double analytical = du_dw(i, j * 2 + i);

                REQUIRE_THAT(analytical, WithinAbs(fd_derivative, 1e-5));
            }
        }
    }
}

TEST_CASE("BarycentricInterpolator - Fixed vs Dynamic Dimensions", "[barycentric][dimensions]")
{
    SECTION("Comparing fixed and dynamic implementations")
    {
        constexpr int N = 6;

        // Create identical nodes
        Eigen::GMatrix<double, N, 1> fixed_nodes;
        Eigen::GMatrix<double, Eigen::Dynamic, 1> dynamic_nodes(N);

        for (int i = 0; i < N; ++i)
        {
            double val = static_cast<double>(i) / (N - 1);
            fixed_nodes(i) = val;
            dynamic_nodes(i) = val;
        }

        // Create interpolators
        BarycentricInterpolatorTpl<double, N> fixed_interp(fixed_nodes);
        BarycentricInterpolatorTpl<double, Eigen::Dynamic> dynamic_interp(dynamic_nodes);

        // Create test values (2×6 matrices)
        Eigen::GMatrix<double, 2, N> fixed_values;
        Eigen::GMatrix<double, 2, Eigen::Dynamic> dynamic_values(2, N);

        for (int i = 0; i < N; ++i)
        {
            double val1 = std::sin(2 * M_PI * fixed_nodes(i));
            double val2 = std::cos(2 * M_PI * fixed_nodes(i));
            fixed_values(0, i) = val1;
            fixed_values(1, i) = val2;
            dynamic_values(0, i) = val1;
            dynamic_values(1, i) = val2;
        }

        // Compare results
        for (double t : {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0})
        {
            Eigen::GMatrix<double, 2, 1> fixed_result, dynamic_result;

            fixed_interp.calc(t, fixed_values, fixed_result);
            dynamic_interp.calc(t, dynamic_values, dynamic_result);

            REQUIRE_THAT(fixed_result(0, 0), WithinAbs(dynamic_result(0, 0), TOLERANCE<double>));
            REQUIRE_THAT(fixed_result(1, 0), WithinAbs(dynamic_result(1, 0), TOLERANCE<double>));
        }
    }

    SECTION("Using DimensionTpl with interpolator")
    {
        // Test with fixed dimension
        constexpr DimensionTpl<4> fixed_dim;
        Eigen::GMatrix<double, fixed_dim.value(), 1> nodes1;
        nodes1 << 0.0, 0.33, 0.67, 1.0;

        BarycentricInterpolatorTpl<double, fixed_dim.Value> interp1(nodes1);
        REQUIRE(interp1.get_n() == 4);
        REQUIRE(interp1.NDim().IsFixed);

        // Test with dynamic dimension
        DimensionTpl<> dynamic_dim(5);
        Eigen::GMatrix<double, Eigen::Dynamic, 1> nodes2(dynamic_dim.value());
        for (int i = 0; i < dynamic_dim.value(); ++i)
        {
            nodes2(i) = static_cast<double>(i) / (dynamic_dim.value() - 1);
        }

        BarycentricInterpolatorTpl<double, Eigen::Dynamic> interp2(nodes2);
        REQUIRE(interp2.get_n() == 5);
        REQUIRE(interp2.NDim().IsDynamic);
    }

    SECTION("Block operations with DimensionTpl")
    {
        constexpr int N = 4;
        constexpr int DIM = 3;

        Eigen::GMatrix<double, N, 1> nodes = Eigen::GMatrix<double, N, 1>::LinSpaced(N, 0.0, 1.0);
        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Create a larger matrix and use block operations
        Eigen::GMatrix<double, DIM + 2, N> full_values;
        full_values.setRandom();

        // Extract a block using DimensionTpl
        constexpr DimensionTpl<DIM> dim_size;
        auto block_values = galileo::block(full_values, 1, 0, dim_size, interp.NDim());

        REQUIRE(block_values.rows() == DIM);
        REQUIRE(block_values.cols() == N);

        // Perform interpolation on the block
        double t = 0.5;
        Eigen::GMatrix<double, DIM, 1> result;
        interp.calc(t, block_values, result);

        REQUIRE(result.rows() == DIM);
        REQUIRE(result.cols() == 1);
    }
}

TEMPLATE_TEST_CASE("BarycentricInterpolator - Different Scalar Types",
                   "[barycentric][scalar_types]",
                   float, double, long double)
{
    using Scalar = TestType;
    constexpr Scalar tolerance = TOLERANCE<Scalar>;
    constexpr int N = 15;

    SECTION("Basic interpolation with different scalar types")
    {
        Eigen::GMatrix<Scalar, N, 1> nodes;
        for (int i = 0; i < N; ++i)
        {
            nodes(i) = static_cast<Scalar>(i) / static_cast<Scalar>(N - 1);
        }

        BarycentricInterpolatorTpl<Scalar, N> interp(nodes);

        // Test with a simple function (1×5 matrix)
        Eigen::GMatrix<Scalar, 1, N> values;
        for (int i = 0; i < N; ++i)
        {
            values(0, i) = nodes(i) * nodes(i); // f(x) = x^2
        }

        // Test interpolation
        for (int i = 0; i <= 10; ++i)
        {
            Scalar t = static_cast<Scalar>(i) / static_cast<Scalar>(10);
            Scalar expected = t * t;

            Eigen::GMatrix<Scalar, 1, 1> result;
            interp.calc(t, values, result);

            REQUIRE_THAT(static_cast<double>(result(0, 0)),
                         WithinAbs(static_cast<double>(expected),
                                   static_cast<double>(tolerance) * 100));
        }
    }

    SECTION("Numerical stability with different precisions")
    {
        // Use Chebyshev nodes for better stability
        Eigen::GMatrix<Scalar, N, 1> nodes;
        for (int i = 0; i < N; ++i)
        {
            Scalar theta = M_PI * (2 * i + 1) / (2 * N);
            nodes(i) = static_cast<Scalar>(0.5) * (Scalar(1.0) - std::cos(theta));
        }

        BarycentricInterpolatorTpl<Scalar, N> interp(nodes);

        // Test with exponential function (2×5 matrix for 2D interpolation)
        Eigen::GMatrix<Scalar, 2, N> values;
        for (int i = 0; i < N; ++i)
        {
            values(0, i) = std::exp(nodes(i));
            values(1, i) = std::sin(static_cast<Scalar>(M_PI) * nodes(i));
        }

        // Compute max error
        Scalar max_error_exp = 0;
        Scalar max_error_sin = 0;
        for (int i = 0; i <= 50; ++i)
        {
            Scalar t = static_cast<Scalar>(i) / static_cast<Scalar>(50);
            Scalar expected_exp = std::exp(t);
            Scalar expected_sin = std::sin(static_cast<Scalar>(M_PI) * t);

            Eigen::GMatrix<Scalar, 2, 1> result;
            interp.calc(t, values, result);

            Scalar error_exp = std::abs(result(0, 0) - expected_exp);
            Scalar error_sin = std::abs(result(1, 0) - expected_sin);
            max_error_exp = std::max(max_error_exp, error_exp);
            max_error_sin = std::max(max_error_sin, error_sin);
        }

        // Expected error depends on precision
        if constexpr (std::is_same_v<Scalar, float>)
        {
            REQUIRE(max_error_exp < Scalar(1e-3));
            REQUIRE(max_error_sin < Scalar(1e-3));
        }
        else if constexpr (std::is_same_v<Scalar, double>)
        {
            REQUIRE(max_error_exp < Scalar(1e-6));
            REQUIRE(max_error_sin < Scalar(1e-6));
        }
        else // long double
        {
            REQUIRE(max_error_exp < Scalar(1e-9));
            REQUIRE(max_error_sin < Scalar(1e-9));
        }
    }
}

TEST_CASE("BarycentricInterpolator - Numerical Stability", "[barycentric][stability]")
{
    SECTION("Clustered nodes")
    {
        constexpr int N = 7;
        Eigen::GMatrix<double, N, 1> nodes;

        // Create clustered nodes near 0
        nodes << 0.0, 0.01, 0.02, 0.03, 0.5, 0.97, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Test with a smooth function (1×7 matrix)
        Eigen::GMatrix<double, 1, N> values;
        for (int i = 0; i < N; ++i)
        {
            values(0, i) = std::sin(M_PI * nodes(i));
        }

        // Should still interpolate reasonably
        for (double t : {0.1, 0.3, 0.5, 0.7, 0.9})
        {
            Eigen::GMatrix<double, 1, 1> result;
            REQUIRE_NOTHROW(interp.calc(t, values, result));

            // Result should be bounded
            REQUIRE(result(0, 0) >= -1.1);
            REQUIRE(result(0, 0) <= 1.1);
        }
    }

    SECTION("Large value ranges")
    {
        constexpr int N = 5;
        Eigen::GMatrix<double, N, 1> nodes = Eigen::GMatrix<double, N, 1>::LinSpaced(N, 0.0, 1.0);

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Test with large values (1×5 matrix)
        Eigen::GMatrix<double, 1, N> values;
        values << 1e6, -2e6, 3e6, -1e6, 2e6;

        // Should handle large values without overflow
        for (double t : {0.2, 0.4, 0.6, 0.8})
        {
            Eigen::GMatrix<double, 1, 1> result;
            REQUIRE_NOTHROW(interp.calc(t, values, result));
            REQUIRE(std::isfinite(result(0, 0)));
        }
    }

    SECTION("Near-singular configurations")
    {
        constexpr int N = 4;
        Eigen::GMatrix<double, N, 1> nodes;

        // Nodes with very small spacing
        nodes << 0.0, 0.5, 0.5 + 1e-10, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        Eigen::GMatrix<double, 1, N> values;
        values << 1.0, 2.0, 2.0 + 1e-8, 3.0;

        // Should still work, though accuracy may be reduced
        double t = 0.7;
        Eigen::GMatrix<double, 1, 1> result;
        REQUIRE_NOTHROW(interp.calc(t, values, result));
        REQUIRE(std::isfinite(result(0, 0)));
    }
}

TEST_CASE("BarycentricInterpolator - Error Handling", "[barycentric][errors]")
{
    SECTION("Invalid t values")
    {
        constexpr int N = 3;
        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.5, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        Eigen::GMatrix<double, 1, N> values;
        values << 1.0, 2.0, 3.0;

        Eigen::GMatrix<double, 1, 1> result;

        // t < 0
        REQUIRE_THROWS_AS(interp.calc(-0.1, values, result), std::runtime_error);

        // t > 1
        REQUIRE_THROWS_AS(interp.calc(1.1, values, result), std::runtime_error);
    }

    SECTION("Dimension mismatch")
    {
        constexpr int N = 4;
        Eigen::GMatrix<double, N, 1> nodes = Eigen::GMatrix<double, N, 1>::LinSpaced(N, 0.0, 1.0);

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Wrong number of columns
        Eigen::GMatrix<double, 2, 3> wrong_values; // Should be 2×4
        wrong_values.setRandom();

        Eigen::GMatrix<double, 2, 1> result;
        REQUIRE_THROWS_AS(interp.calc(0.5, wrong_values, result), std::runtime_error);
    }

    SECTION("calcDiff dimension requirements")
    {
        constexpr int N = 3;
        Eigen::GMatrix<double, N, 1> nodes;
        nodes << 0.0, 0.5, 1.0;

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        Eigen::GMatrix<double, 2, N> values; // 2×3 matrix
        values.setRandom();

        // Wrong output dimensions
        Eigen::GMatrix<double, 2, 5> wrong_du_dw; // Should be 2×6
        REQUIRE_THROWS_AS(interp.calcDiff(0.5, values, wrong_du_dw), std::runtime_error);

        // Correct dimensions
        Eigen::GMatrix<double, 2, 6> correct_du_dw; // 2×(2*3)
        REQUIRE_NOTHROW(interp.calcDiff(0.5, values, correct_du_dw));
    }

    SECTION("Fixed dimension mismatch in constructor")
    {
        // Try to construct with wrong number of nodes
        Eigen::GMatrix<double, 5, 1> nodes = Eigen::GMatrix<double, 5, 1>::LinSpaced(5, 0.0, 1.0);

        // This should fail at runtime due to dimension mismatch
        REQUIRE_THROWS_AS((BarycentricInterpolatorTpl<double, 3>(nodes)), std::runtime_error);
    }
}

TEST_CASE("BarycentricInterpolator - Performance", "[barycentric][performance]")
{
    SECTION("Large interpolation problem")
    {
        constexpr int N = 50;
        constexpr int M = 100;     // Number of output dimensions
        constexpr int NPTS = 1000; // Number of evaluation points

        // Use Chebyshev nodes for good conditioning
        Eigen::GMatrix<double, N, 1> nodes;
        for (int i = 0; i < N; ++i)
        {
            double theta = M_PI * (2 * i + 1) / (2 * N);
            nodes(i) = 0.5 * (1.0 - std::cos(theta));
        }

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Create large data matrix (100×50)
        Eigen::GMatrix<double, M, N> values;
        values.setRandom();

        // Time interpolation
        auto start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < NPTS; ++i)
        {
            double t = static_cast<double>(i) / (NPTS - 1);
            Eigen::GMatrix<double, M, 1> result;
            interp.calc(t, values, result);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        // Should complete in reasonable time
        REQUIRE(duration.count() < 1000); // Less than 1 second for 1000 evaluations
    }

    SECTION("Comparison of fixed vs dynamic performance")
    {
        constexpr int N = 10;
        constexpr int NITER = 10000;

        // Fixed dimension setup
        Eigen::GMatrix<double, N, 1> fixed_nodes = Eigen::GMatrix<double, N, 1>::LinSpaced(N, 0.0, 1.0);
        BarycentricInterpolatorTpl<double, N> fixed_interp(fixed_nodes);

        // Dynamic dimension setup
        Eigen::GMatrix<double, Eigen::Dynamic, 1> dynamic_nodes = Eigen::GMatrix<double, Eigen::Dynamic, 1>::LinSpaced(N, 0.0, 1.0);
        BarycentricInterpolatorTpl<double, Eigen::Dynamic> dynamic_interp(dynamic_nodes);

        Eigen::GMatrix<double, 3, N> fixed_values;
        Eigen::GMatrix<double, 3, Eigen::Dynamic> dynamic_values(3, N);
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                double val = std::sin(2 * M_PI * i / (N - 1) + j);
                fixed_values(j, i) = val;
                dynamic_values(j, i) = val;
            }
        }

        // Time fixed dimension
        auto start_fixed = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < NITER; ++i)
        {
            double t = 0.5;
            Eigen::GMatrix<double, 3, 1> result;
            fixed_interp.calc(t, fixed_values, result);
        }
        auto end_fixed = std::chrono::high_resolution_clock::now();

        // Time dynamic dimension
        auto start_dynamic = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < NITER; ++i)
        {
            double t = 0.5;
            Eigen::GMatrix<double, 3, 1> result;
            dynamic_interp.calc(t, dynamic_values, result);
        }
        auto end_dynamic = std::chrono::high_resolution_clock::now();

        auto fixed_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_fixed - start_fixed);
        auto dynamic_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_dynamic - start_dynamic);

        // Fixed dimension should generally be faster due to compile-time optimizations
        // But we just check that both complete in reasonable time
        REQUIRE(fixed_duration.count() > 0);
        REQUIRE(dynamic_duration.count() > 0);
    }
}

TEMPLATE_TEST_CASE("BarycentricInterpolator - Integration with JacobiRoots, different nodes", "[barycentric][Jacobi params]",
                   LegendreParams, Chebyshev1Params, Chebyshev2Params)
{
    SECTION("Optimal interpolation with Jacobi roots")
    {
        constexpr int N = 100;

        // Different Jacobi parameters

        JacobiRootsTpl<double, N> jacobi_roots(TestType::ALPHA, TestType::BETA);
        jacobi_roots.compute_roots();
        const auto &nodes = jacobi_roots.get_roots();

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Test with Runge function (1×N matrix)
        auto runge = [](double x)
        {
            double shifted = 2.0 * x - 1.0;
            return 1.0 / (1.0 + 25.0 * shifted * shifted);
        };

        Eigen::GMatrix<double, 1, N> values;
        for (int i = 0; i < N; ++i)
        {
            values(0, i) = runge(nodes(i));
        }

        // Compute max error
        double max_error = 0.0;
        int n_points = 1000;
        for (int i = 0; i <= n_points; ++i)
        {
            double t = i / static_cast<double>(n_points);
            double expected = runge(t);

            Eigen::GMatrix<double, 1, 1> result;
            interp.calc(t, values, result);

            max_error = std::max(max_error, std::abs(result(0, 0) - expected));
        }

        // With enough nodes, even the notoriously challenging Runge function should be manageable
        REQUIRE(max_error < 1e-7);
    }

    SECTION("High-order interpolation with Jacobi roots")
    {
        constexpr int N = 20;

        // Use Legendre nodes for high-order interpolation
        JacobiRootsTpl<double, N> jacobi_roots(0.0, 0.0);
        jacobi_roots.compute_roots();
        const auto &nodes = jacobi_roots.get_roots();

        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Test with smooth function (2×15 matrix for 2D)
        auto func1 = [](double x)
        {
            return std::exp(x) * std::sin(4 * M_PI * x);
        };
        auto func2 = [](double x)
        {
            return std::cos(3 * M_PI * x) / (1.0 + x);
        };

        Eigen::GMatrix<double, 2, N> values;
        for (int i = 0; i < N; ++i)
        {
            values(0, i) = func1(nodes(i));
            values(1, i) = func2(nodes(i));
        }

        // Should achieve high accuracy for smooth functions
        double max_error1 = 0.0;
        double max_error2 = 0.0;
        for (int i = 0; i <= 200; ++i)
        {
            double t = i / 200.0;
            double expected1 = func1(t);
            double expected2 = func2(t);

            Eigen::GMatrix<double, 2, 1> result;
            interp.calc(t, values, result);

            max_error1 = std::max(max_error1, std::abs(result(0, 0) - expected1));
            max_error2 = std::max(max_error2, std::abs(result(1, 0) - expected2));
        }

        REQUIRE(max_error1 < 1e-6); // High accuracy for smooth functions
        REQUIRE(max_error2 < 1e-6);
    }
}

TEST_CASE("BarycentricInterpolator - Debug Validation", "[barycentric][debug]")
{
    SECTION("Simple linear interpolation validation")
    {
        // Test with 2 nodes - should be exact linear interpolation
        Eigen::GMatrix<double, 2, 1> nodes;
        nodes << 0.0, 1.0;

        BarycentricInterpolatorTpl<double, 2> interp(nodes);

        // Linear function f(x) = 2x + 1
        Eigen::GMatrix<double, 1, 2> values;
        values << 1.0, 3.0; // f(0) = 1, f(1) = 3

        // Test at various points - should be exact
        std::vector<double> test_points = {0.0, 0.25, 0.5, 0.75, 1.0};
        for (double t : test_points)
        {
            Eigen::GMatrix<double, 1, 1> result;
            interp.calc(t, values, result);
            double expected = 1.0 + 2.0 * t;

            REQUIRE_THAT(result(0, 0), WithinAbs(expected, 1e-14));
        }
    }

    SECTION("Debug vector interpolation")
    {
        // Use the same setup as failing test but with simpler function
        constexpr int N = 5;
        constexpr int DIM = 3;

        Eigen::GMatrix<double, N, 1> nodes = Eigen::GMatrix<double, N, 1>::LinSpaced(N, 0.0, 1.0);
        BarycentricInterpolatorTpl<double, N> interp(nodes);

        // Create simple linear data instead of trigonometric
        Eigen::GMatrix<double, DIM, N> values;
        for (int i = 0; i < N; ++i)
        {
            values(0, i) = 2.0 * nodes(i) + 1.0; // f(x) = 2x + 1
            values(1, i) = 3.0 * nodes(i) - 0.5; // f(x) = 3x - 0.5
            values(2, i) = nodes(i) * nodes(i);  // f(x) = x²
        }

        // Test at t = 0.1
        double t = 0.1;
        Eigen::GMatrix<double, DIM, 1> result;
        interp.calc(t, values, result);

        double expected0 = 2.0 * t + 1.0;
        double expected1 = 3.0 * t - 0.5;
        double expected2 = t * t;

        REQUIRE_THAT(result(0, 0), WithinAbs(expected0, 1e-14));
        REQUIRE_THAT(result(1, 0), WithinAbs(expected1, 1e-14));
        REQUIRE_THAT(result(2, 0), WithinAbs(expected2, 1e-14));
    }
}
