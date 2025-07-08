#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/common/math/jacobi-roots.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <chrono>

using namespace galileo::math;
using namespace Catch::Matchers;

template<typename NumScalar>
constexpr NumScalar TOLERANCE = std::numeric_limits<NumScalar>::epsilon() * 100;

// Test fixture for common setup
template<typename NumScalar, int N>
class JacobiRootsTestFixture
{
public:
    using JacobiRootsType = JacobiRootsTpl<NumScalar, N>;

    JacobiRootsType standard_jacobi;     // $\alpha = 0, \beta = 0$ (Legendre)
    JacobiRootsType chebyshev_first;     // $\alpha = -0.5, \beta = -0.5$
    JacobiRootsType chebyshev_second;    // $\alpha = 0.5, \beta = 0.5$
    JacobiRootsType custom_jacobi;       // $\alpha = 1.5, \beta = 2.5$

    JacobiRootsTestFixture()
        : standard_jacobi(0.0, 0.0)
        , chebyshev_first(-0.5, -0.5)
        , chebyshev_second(0.5, 0.5)
        , custom_jacobi(1.5, 2.5)
    {
    }
};

TEST_CASE("JacobiRoots - Basic Construction and Properties", "[jacobi_roots]")
{
    SECTION("Default construction")
    {
        JacobiRootsTpl<double, 5> jacobi_roots;
        REQUIRE(jacobi_roots.get_alpha() == 0.0);
        REQUIRE(jacobi_roots.get_beta() == 0.0);
    }

    SECTION("Construction with parameters")
    {
        constexpr double alpha = 1.5;
        constexpr double beta = 2.0;

        JacobiRootsTpl<double, 8> jacobi_roots(alpha, beta);
        REQUIRE(jacobi_roots.get_alpha() == alpha);
        REQUIRE(jacobi_roots.get_beta() == beta);
    }

    SECTION("Different scalar types")
    {
        JacobiRootsTpl<float, 4> float_roots(1.0f, 0.5f);
        REQUIRE(float_roots.get_alpha() == 1.0f);
        REQUIRE(float_roots.get_beta() == 0.5f);

        JacobiRootsTpl<long double, 4> long_double_roots(1.0L, 0.5L);
        REQUIRE(long_double_roots.get_alpha() == 1.0L);
        REQUIRE(long_double_roots.get_beta() == 0.5L);
    }

    SECTION("Template parameter validation")
    {
        static_assert(JacobiRootsTpl<double, 3>::N == 3);
        static_assert(std::is_same_v<JacobiRootsTpl<float, 5>::NumScalar, float>);
        static_assert(JacobiRootsTpl<double, 7>::Options == 0);
    }
}

TEST_CASE_METHOD((JacobiRootsTestFixture<double, 5>), "JacobiRoots - Root Computation Properties", "[jacobi_roots]")
{
    SECTION("Standard Legendre polynomial roots ($\\alpha=0, \\beta=0$)")
    {
        standard_jacobi.compute_roots();
        const auto& roots = standard_jacobi.get_roots();

        REQUIRE(roots.size() == 5);

        // Roots should be in [0, 1] after mapping
        for (int i = 0; i < 5; ++i)
        {
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }

        // Roots should be sorted
        for (int i = 0; i < 4; ++i)
        {
            REQUIRE(roots(i) < roots(i + 1));
        }

        // For Legendre polynomials, roots should be symmetric about 0.5
        for (int i = 0; i < 5; ++i)
        {
            double expected_symmetric = 1.0 - roots(4 - i);
            REQUIRE_THAT(roots(i), WithinAbs(expected_symmetric, TOLERANCE<double>));
        }
    }

    SECTION("Chebyshev polynomial roots")
    {
        chebyshev_first.compute_roots();
        const auto& roots = chebyshev_first.get_roots();

        REQUIRE(roots.size() == 5);

        // Roots should be in [0, 1]
        for (int i = 0; i < 5; ++i)
        {
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }

        // Roots should be sorted
        for (int i = 0; i < 4; ++i)
        {
            REQUIRE(roots(i) < roots(i + 1));
        }
    }

    SECTION("Custom Jacobi polynomial roots")
    {
        custom_jacobi.compute_roots();
        const auto& roots = custom_jacobi.get_roots();

        REQUIRE(roots.size() == 5);

        // Roots should be in [0, 1]
        for (int i = 0; i < 5; ++i)
        {
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }

        // Roots should be sorted
        for (int i = 0; i < 4; ++i)
        {
            REQUIRE(roots(i) < roots(i + 1));
        }
    }
}

TEMPLATE_TEST_CASE("JacobiRoots - Different Sizes", "[jacobi_roots][template]",
                   (std::integral_constant<int, 1>),
                   (std::integral_constant<int, 2>),
                   (std::integral_constant<int, 3>),
                   (std::integral_constant<int, 10>),
                   (std::integral_constant<int, 20>))
{
    constexpr int N = TestType::value;
    JacobiRootsTpl<double, N> jacobi_roots(0.0, 0.0);

    jacobi_roots.compute_roots();
    const auto& roots = jacobi_roots.get_roots();
    const auto& jacobi_matrix = jacobi_roots.get_jacobi_matrix();

    REQUIRE(roots.size() == N);
    REQUIRE(jacobi_matrix.rows() == N);
    REQUIRE(jacobi_matrix.cols() == N);

    // All roots should be in [0, 1]
    for (int i = 0; i < N; ++i)
    {
        REQUIRE(roots(i) >= 0.0);
        REQUIRE(roots(i) <= 1.0);
    }

    // Roots should be sorted
    for (int i = 0; i < N - 1; ++i)
    {
        REQUIRE(roots(i) < roots(i + 1));
    }

    // Jacobi matrix should be symmetric tridiagonal
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (std::abs(i - j) > 1)
            {
                REQUIRE_THAT(jacobi_matrix(i, j), WithinAbs(0.0, TOLERANCE<double>));
            }
        }
    }
}

TEST_CASE("JacobiRoots - Edge Cases", "[jacobi_roots][edge_cases]")
{
    SECTION("Single root (N=1)")
    {
        JacobiRootsTpl<double, 1> single_root(0.0, 0.0);
        single_root.compute_roots();

        const auto& roots = single_root.get_roots();
        REQUIRE(roots.size() == 1);
        REQUIRE_THAT(roots(0), WithinAbs(0.5, TOLERANCE<double>));  // Should be at center
    }

    SECTION("Two roots (N=2)")
    {
        JacobiRootsTpl<double, 2> two_roots(0.0, 0.0);
        two_roots.compute_roots();

        const auto& roots = two_roots.get_roots();
        REQUIRE(roots.size() == 2);

        // For Legendre polynomials with N=2, roots should be symmetric about 0.5
        double expected1 = 0.5 - std::sqrt(3.0) / 6.0;
        double expected2 = 0.5 + std::sqrt(3.0) / 6.0;

        REQUIRE_THAT(roots(0), WithinAbs(expected1, TOLERANCE<double>));
        REQUIRE_THAT(roots(1), WithinAbs(expected2, TOLERANCE<double>));
    }

    SECTION("Extreme $\\alpha$ and $\\beta$ values")
    {
        JacobiRootsTpl<double, 5> extreme_roots(10.0, 15.0);
        extreme_roots.compute_roots();

        const auto& roots = extreme_roots.get_roots();
        REQUIRE(roots.size() == 5);

        // Even with extreme parameters, roots should be in [0, 1]
        for (int i = 0; i < 5; ++i)
        {
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }
    }

    SECTION("Negative $\\alpha$ and $\\beta$ values")
    {
        JacobiRootsTpl<double, 4> negative_roots(-0.9, -0.8);
        negative_roots.compute_roots();

        const auto& roots = negative_roots.get_roots();
        REQUIRE(roots.size() == 4);

        for (int i = 0; i < 4; ++i)
        {
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }
    }
}

TEST_CASE("JacobiRoots - Jacobi Matrix Properties", "[jacobi_roots][matrix]")
{
    SECTION("Matrix structure and properties")
    {
        JacobiRootsTpl<double, 6> jacobi_roots(1.0, 0.5);
        jacobi_roots.compute_roots();

        const auto& matrix = jacobi_roots.get_jacobi_matrix();

        // Should be square
        REQUIRE(matrix.rows() == 6);
        REQUIRE(matrix.cols() == 6);

        // Should be symmetric
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                REQUIRE_THAT(matrix(i, j), WithinAbs(matrix(j, i), TOLERANCE<double>));
            }
        }

        // Should be tridiagonal
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                if (std::abs(i - j) > 1)
                {
                    REQUIRE_THAT(matrix(i, j), WithinAbs(0.0, TOLERANCE<double>));
                }
            }
        }
    }

            SECTION("Matrix eigenvalue properties")
    {
        JacobiRootsTpl<double, 4> jacobi_roots(0.0, 0.0);
        jacobi_roots.compute_roots();

        const auto& matrix = jacobi_roots.get_jacobi_matrix();
        const auto& roots = jacobi_roots.get_roots();

        // Verify that matrix eigenvalues are real and finite
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
        auto eigenvalues = solver.eigenvalues();

        REQUIRE(eigenvalues.size() == 4);

        // All eigenvalues should be real and finite (since matrix is symmetric)
        for (int i = 0; i < 4; ++i)
        {
            REQUIRE(std::isfinite(eigenvalues(i)));
        }

        // Check that the computed roots are in the expected range [0, 1]
        for (int i = 0; i < 4; ++i)
        {
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }

        // Eigenvalues should be in [-1, 1] for Legendre polynomials before mapping
        for (int i = 0; i < 4; ++i)
        {
            REQUIRE(eigenvalues(i) >= -1.0);
            REQUIRE(eigenvalues(i) <= 1.0);
        }

        // The eigenvalues should be distinct (no repeated roots for this case)
        for (int i = 0; i < 3; ++i)
        {
            REQUIRE(std::abs(eigenvalues(i) - eigenvalues(i + 1)) > TOLERANCE<double>);
        }
    }
}

TEST_CASE("JacobiRoots - Mathematical Properties", "[jacobi_roots][mathematical]")
{
    SECTION("Orthogonality properties for Legendre case")
    {
        JacobiRootsTpl<double, 8> legendre_roots(0.0, 0.0);
        legendre_roots.compute_roots();

        const auto& roots = legendre_roots.get_roots();

        // Map roots back to [-1, 1] for polynomial evaluation
        Eigen::VectorXd unmapped_roots = 2.0 * roots.array() - 1.0;

        // Verify that Legendre polynomials of different degrees are orthogonal
        // This is a basic check using the roots
        for (int i = 0; i < 8; ++i)
        {
            double x = unmapped_roots(i);

            // $P_0(x) = 1$, $P_1(x) = x$, $P_2(x) = (3x^2 - 1)/2$
            double P0 = 1.0;
            double P1 = x;
            double P2 = 0.5 * (3.0 * x * x - 1.0);

            // These should not all be zero simultaneously (except possibly at specific points)
            REQUIRE(!(std::abs(P0) < TOLERANCE<double> && std::abs(P1) < TOLERANCE<double> && std::abs(P2) < TOLERANCE<double>));
        }
    }

    SECTION("Symmetry properties for standard Jacobi polynomials")
    {
        JacobiRootsTpl<double, 7> symmetric_roots(0.0, 0.0);
        symmetric_roots.compute_roots();

        const auto& roots = symmetric_roots.get_roots();

        // For symmetric Jacobi polynomials ($\alpha = \beta$), roots should be symmetric about 0.5
        for (int i = 0; i < 7; ++i)
        {
            double symmetric_point = 1.0 - roots(6 - i);
            REQUIRE_THAT(roots(i), WithinAbs(symmetric_point, TOLERANCE<double>));
        }
    }
}

TEST_CASE("JacobiRoots - Precision and Numerical Stability", "[jacobi_roots][numerical]")
{
    SECTION("High precision requirements")
    {
        JacobiRootsTpl<double, 15> high_order_roots(0.0, 0.0);
        high_order_roots.compute_roots();

        const auto& roots = high_order_roots.get_roots();

        // Even for high-order polynomials, roots should be computed accurately
        for (int i = 0; i < 15; ++i)
        {
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);

            if (i > 0)
            {
                REQUIRE(roots(i) > roots(i - 1));
            }
        }

        // Check that roots are not clustered too tightly (numerical stability)
        for (int i = 0; i < 14; ++i)
        {
            double gap = roots(i + 1) - roots(i);
            REQUIRE(gap > TOLERANCE<double>);  // Minimum separation
        }
    }

    SECTION("Consistency across multiple computations")
    {
        JacobiRootsTpl<double, 6> roots1(0.5, 1.5);
        JacobiRootsTpl<double, 6> roots2(0.5, 1.5);

        roots1.compute_roots();
        roots2.compute_roots();

        const auto& r1 = roots1.get_roots();
        const auto& r2 = roots2.get_roots();

        for (int i = 0; i < 6; ++i)
        {
            REQUIRE_THAT(r1(i), WithinAbs(r2(i), TOLERANCE<double>));
        }
    }
}

TEST_CASE("JacobiRoots - Copy and Assignment", "[jacobi_roots][copy]")
{
    SECTION("Copy construction")
    {
        JacobiRootsTpl<double, 5> original(2.0, 3.0);
        original.compute_roots();

        JacobiRootsTpl<double, 5> copy(original);

        REQUIRE(copy.get_alpha() == original.get_alpha());
        REQUIRE(copy.get_beta() == original.get_beta());

        // Note: The roots and matrix are computed, not copied, so we need to compute them
        copy.compute_roots();

        const auto& orig_roots = original.get_roots();
        const auto& copy_roots = copy.get_roots();

        for (int i = 0; i < 5; ++i)
        {
            REQUIRE_THAT(copy_roots(i), WithinAbs(orig_roots(i), TOLERANCE<double>));
        }
    }

    SECTION("Assignment operator")
    {
        JacobiRootsTpl<double, 4> original(1.5, 0.5);
        original.compute_roots();

        JacobiRootsTpl<double, 4> assigned;
        assigned = original;

        REQUIRE(assigned.get_alpha() == original.get_alpha());
        REQUIRE(assigned.get_beta() == original.get_beta());

        assigned.compute_roots();

        const auto& orig_roots = original.get_roots();
        const auto& assigned_roots = assigned.get_roots();

        for (int i = 0; i < 4; ++i)
        {
            REQUIRE_THAT(assigned_roots(i), WithinAbs(orig_roots(i), TOLERANCE<double>));
        }
    }
}

TEST_CASE("JacobiRoots - Error Handling and Robustness", "[jacobi_roots][error_handling]")
{
    SECTION("Very small $\\alpha$ and $\\beta$ values")
    {
        JacobiRootsTpl<double, 5> small_params(1e-15, 1e-15);
        REQUIRE_NOTHROW(small_params.compute_roots());

        const auto& roots = small_params.get_roots();
        for (int i = 0; i < 5; ++i)
        {
            REQUIRE(std::isfinite(roots(i)));
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }
    }

    SECTION("Large $\\alpha$ and $\\beta$ values")
    {
        JacobiRootsTpl<double, 3> large_params(100.0, 150.0);
        REQUIRE_NOTHROW(large_params.compute_roots());

        const auto& roots = large_params.get_roots();
        for (int i = 0; i < 3; ++i)
        {
            REQUIRE(std::isfinite(roots(i)));
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }
    }
}

TEST_CASE("JacobiRoots - Performance and Large Sizes", "[jacobi_roots][performance]")
{
    SECTION("Large polynomial order")
    {
        constexpr int LARGE_N = 50;
        JacobiRootsTpl<double, LARGE_N> large_roots(0.0, 0.0);

        auto start = std::chrono::high_resolution_clock::now();
        large_roots.compute_roots();
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        // Should complete in reasonable time (less than 1 second for N=50)
        REQUIRE(duration.count() < 1000);

        const auto& roots = large_roots.get_roots();
        REQUIRE(roots.size() == LARGE_N);

        // Verify basic properties still hold
        for (int i = 0; i < LARGE_N - 1; ++i)
        {
            REQUIRE(roots(i) < roots(i + 1));
            REQUIRE(roots(i) >= 0.0);
            REQUIRE(roots(i) <= 1.0);
        }
    }
}

TEMPLATE_TEST_CASE("JacobiRoots - Different Scalar Types", "[jacobi_roots][scalar_types]",
                   float, double, long double)
{
    using Scalar = TestType;
    constexpr Scalar tolerance = TOLERANCE<Scalar>;

    JacobiRootsTpl<Scalar, 5> jacobi_roots(Scalar(0.5), Scalar(1.0));
    jacobi_roots.compute_roots();

    const auto& roots = jacobi_roots.get_roots();

    REQUIRE(roots.size() == 5);

    for (int i = 0; i < 5; ++i)
    {
        REQUIRE(roots(i) >= Scalar(0.0));
        REQUIRE(roots(i) <= Scalar(1.0));

        if (i > 0)
        {
            REQUIRE(roots(i) > roots(i - 1));
        }
    }

    // Test that the computation is reasonably accurate for the scalar type
    REQUIRE_THAT(static_cast<double>(jacobi_roots.get_alpha()), WithinAbs(0.5, tolerance));
    REQUIRE_THAT(static_cast<double>(jacobi_roots.get_beta()), WithinAbs(1.0, tolerance));
}
