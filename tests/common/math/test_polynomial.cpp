#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/common/math/polynomial.hpp"
#include <Eigen/Dense>
#include <cmath>

using namespace galileo::math;
using namespace Catch::Matchers;

// Test fixture for common polynomial setups
template<typename Scalar, int N>
class PolynomialTestFixture
{
public:
    static constexpr int PolynomialDegree = N;
    using PolynomialType = JacobiPolynomialTpl<Scalar, N, Eigen::ColMajor>;
    using ScalarType = Scalar;

    // Common alpha, beta values for testing
    static constexpr Scalar alpha_legendre = Scalar(0.0);
    static constexpr Scalar beta_legendre = Scalar(0.0);
    static constexpr Scalar alpha_chebyshev_first = Scalar(-0.5);
    static constexpr Scalar beta_chebyshev_first = Scalar(-0.5);
    static constexpr Scalar alpha_chebyshev_second = Scalar(0.5);
    static constexpr Scalar beta_chebyshev_second = Scalar(0.5);

    PolynomialTestFixture() = default;
};

TEST_CASE("JacobiPolynomial - Basic Construction and Properties", "[polynomial]")
{
    SECTION("Default construction")
    {
        JacobiPolynomialTpl<double, 5> poly;
        // Default constructor should set alpha = 0, beta = 0
        // Note: The default constructor doesn't compute nodes/weights
    }

    SECTION("Construction with Legendre parameters (alpha=0, beta=0)")
    {
        constexpr double alpha = 0.0;
        constexpr double beta = 0.0;
        JacobiPolynomialTpl<double, 5> legendre_poly(alpha, beta);

        const auto& nodes = legendre_poly.get_nodes();
        REQUIRE(nodes.size() == 5);

        // All nodes should be in [0, 1] range
        for (int i = 0; i < nodes.size(); ++i)
        {
            REQUIRE(nodes[i] >= 0.0);
            REQUIRE(nodes[i] <= 1.0);
        }

        // Nodes should be sorted in ascending order
        for (int i = 0; i < nodes.size() - 1; ++i)
        {
            REQUIRE(nodes[i] < nodes[i + 1]);
        }
    }

    SECTION("Construction with Chebyshev First Kind parameters (alpha=-0.5, beta=-0.5)")
    {
        constexpr double alpha = -0.5;
        constexpr double beta = -0.5;
        JacobiPolynomialTpl<double, 4> chebyshev_poly(alpha, beta);

        const auto& nodes = chebyshev_poly.get_nodes();
        REQUIRE(nodes.size() == 4);

        // All nodes should be in [0, 1] range
        for (int i = 0; i < nodes.size(); ++i)
        {
            REQUIRE(nodes[i] >= 0.0);
            REQUIRE(nodes[i] <= 1.0);
        }
    }

    SECTION("Construction with Chebyshev Second Kind parameters (alpha=0.5, beta=0.5)")
    {
        constexpr double alpha = 0.5;
        constexpr double beta = 0.5;
        JacobiPolynomialTpl<double, 6> chebyshev2_poly(alpha, beta);

        const auto& nodes = chebyshev2_poly.get_nodes();
        REQUIRE(nodes.size() == 6);

        // All nodes should be in [0, 1] range
        for (int i = 0; i < nodes.size(); ++i)
        {
            REQUIRE(nodes[i] >= 0.0);
            REQUIRE(nodes[i] <= 1.0);
        }
    }

    SECTION("Construction with custom parameters")
    {
        constexpr double alpha = 1.5;
        constexpr double beta = 2.3;
        JacobiPolynomialTpl<double, 3> custom_poly(alpha, beta);

        const auto& nodes = custom_poly.get_nodes();
        REQUIRE(nodes.size() == 3);

        // All nodes should be in [0, 1] range
        for (int i = 0; i < nodes.size(); ++i)
        {
            REQUIRE(nodes[i] >= 0.0);
            REQUIRE(nodes[i] <= 1.0);
        }
    }
}

TEST_CASE("JacobiPolynomial - Node Properties", "[polynomial]")
{
    SECTION("Single node polynomial (N=1)")
    {
        JacobiPolynomialTpl<double, 1> single_poly(0.0, 0.0);
        const auto& nodes = single_poly.get_nodes();

        REQUIRE(nodes.size() == 1);
        REQUIRE_THAT(nodes[0], WithinAbs(0.5, 1e-10)); // Single node should be at center
    }

    SECTION("Two node polynomial (N=2)")
    {
        JacobiPolynomialTpl<double, 2> double_poly(0.0, 0.0);
        const auto& nodes = double_poly.get_nodes();

        REQUIRE(nodes.size() == 2);
        REQUIRE(nodes[0] < 0.5);
        REQUIRE(nodes[1] > 0.5);
        // For Legendre, nodes should be symmetric about 0.5
        REQUIRE_THAT(nodes[0] + nodes[1], WithinAbs(1.0, 1e-10));
    }

    SECTION("Higher order polynomial node count")
    {
        JacobiPolynomialTpl<double, 10> high_order_poly(0.0, 0.0);
        const auto& nodes = high_order_poly.get_nodes();

        REQUIRE(nodes.size() == 10);

        // Check that all nodes are distinct
        for (int i = 0; i < nodes.size(); ++i)
        {
            for (int j = i + 1; j < nodes.size(); ++j)
            {
                REQUIRE(std::abs(nodes[i] - nodes[j]) > 1e-12);
            }
        }
    }

    SECTION("Node stability with different scalar types")
    {
        JacobiPolynomialTpl<float, 5> float_poly(0.0f, 0.0f);
        JacobiPolynomialTpl<double, 5> double_poly(0.0, 0.0);

        const auto& float_nodes = float_poly.get_nodes();
        const auto& double_nodes = double_poly.get_nodes();

        REQUIRE(float_nodes.size() == double_nodes.size());

        // Nodes should be approximately the same (within float precision)
        for (int i = 0; i < float_nodes.size(); ++i)
        {
            REQUIRE_THAT(static_cast<double>(float_nodes[i]), WithinAbs(double_nodes[i], 1e-6));
        }
    }
}

TEST_CASE("JacobiPolynomial - Barycentric Interpolation", "[polynomial]")
{
    SECTION("Interpolation at nodes should return exact values")
    {
        JacobiPolynomialTpl<double, 4> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Create test data: function values at nodes
        Eigen::Matrix<double, 3, 4> w;
        w << 1.0, 2.0, 3.0, 4.0,    // First function
             0.5, 1.5, 2.5, 3.5,    // Second function
             2.0, 4.0, 6.0, 8.0;    // Third function

        Eigen::Vector3d result;

        // Test interpolation at each node
        for (int i = 0; i < nodes.size(); ++i)
        {
            poly.barycentricInterpolation(nodes[i], w, result);

            REQUIRE_THAT(result[0], WithinAbs(w(0, i), 1e-12));
            REQUIRE_THAT(result[1], WithinAbs(w(1, i), 1e-12));
            REQUIRE_THAT(result[2], WithinAbs(w(2, i), 1e-12));
        }
    }

    SECTION("Interpolation of constant function")
    {
        JacobiPolynomialTpl<double, 5> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Constant function: all values are the same
        Eigen::Matrix<double, 2, 5> w;
        w.row(0).setConstant(3.14);
        w.row(1).setConstant(-2.71);

        Eigen::Vector2d result;

        // Test at various points
        std::vector<double> test_points = {0.0, 0.25, 0.5, 0.75, 1.0};
        for (double t : test_points)
        {
            poly.barycentricInterpolation(t, w, result);

            REQUIRE_THAT(result[0], WithinAbs(3.14, 1e-10));
            REQUIRE_THAT(result[1], WithinAbs(-2.71, 1e-10));
        }
    }

    SECTION("Interpolation of linear function")
    {
        JacobiPolynomialTpl<double, 3> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Linear function: f(x) = 2*x + 1
        Eigen::Matrix<double, 1, 3> w;
        for (int i = 0; i < 3; ++i)
        {
            w(0, i) = 2.0 * nodes[i] + 1.0;
        }

        Eigen::Vector1d result;

        // Test at various points - should be exact for polynomials of degree <= N-1
        std::vector<double> test_points = {0.1, 0.3, 0.7, 0.9};
        for (double t : test_points)
        {
            poly.barycentricInterpolation(t, w, result);
            double expected = 2.0 * t + 1.0;

            REQUIRE_THAT(result[0], WithinAbs(expected, 1e-10));
        }
    }

    SECTION("Interpolation of quadratic function")
    {
        JacobiPolynomialTpl<double, 4> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Quadratic function: f(x) = x^2 - 2*x + 3
        Eigen::Matrix<double, 1, 4> w;
        for (int i = 0; i < 4; ++i)
        {
            double x = nodes[i];
            w(0, i) = x * x - 2.0 * x + 3.0;
        }

        Eigen::Vector1d result;

        // Test at various points - should be exact for polynomials of degree <= N-1
        std::vector<double> test_points = {0.15, 0.35, 0.65, 0.85};
        for (double t : test_points)
        {
            poly.barycentricInterpolation(t, w, result);
            double expected = t * t - 2.0 * t + 3.0;

            REQUIRE_THAT(result[0], WithinAbs(expected, 1e-10));
        }
    }

    SECTION("Boundary value interpolation")
    {
        JacobiPolynomialTpl<double, 4> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Simple test function
        Eigen::Matrix<double, 2, 4> w;
        for (int i = 0; i < 4; ++i)
        {
            w(0, i) = std::sin(nodes[i] * M_PI);
            w(1, i) = std::cos(nodes[i] * M_PI);
        }

        Eigen::Vector2d result;

        // Test at boundaries
        poly.barycentricInterpolation(0.0, w, result);
        poly.barycentricInterpolation(1.0, w, result);

        // Should not throw and should produce finite results
        REQUIRE(std::isfinite(result[0]));
        REQUIRE(std::isfinite(result[1]));
    }
}

TEST_CASE("JacobiPolynomial - Barycentric Interpolation Differentiation", "[polynomial]")
{
    SECTION("Differentiation at nodes gives identity for matching columns")
    {
        JacobiPolynomialTpl<double, 3> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Create test data
        Eigen::Matrix<double, 2, 3> w;
        w << 1.0, 2.0, 3.0,
             4.0, 5.0, 6.0;

        Eigen::Matrix<double, 2, 6> du_dw; // 2 x (2*3)

        // Test differentiation at each node
        for (int node_idx = 0; node_idx < nodes.size(); ++node_idx)
        {
            poly.barycentricInterpolationDiff(nodes[node_idx], w, du_dw);

            // Check that the derivative matrix has the right structure
            for (int col = 0; col < 3; ++col)
            {
                if (col == node_idx)
                {
                    // Identity block for the matching node
                    REQUIRE_THAT(du_dw(0, col * 2 + 0), WithinAbs(1.0, 1e-12));
                    REQUIRE_THAT(du_dw(1, col * 2 + 1), WithinAbs(1.0, 1e-12));
                    REQUIRE_THAT(du_dw(0, col * 2 + 1), WithinAbs(0.0, 1e-12));
                    REQUIRE_THAT(du_dw(1, col * 2 + 0), WithinAbs(0.0, 1e-12));
                }
                else
                {
                    // Zero block for non-matching nodes
                    REQUIRE_THAT(du_dw(0, col * 2 + 0), WithinAbs(0.0, 1e-12));
                    REQUIRE_THAT(du_dw(1, col * 2 + 1), WithinAbs(0.0, 1e-12));
                    REQUIRE_THAT(du_dw(0, col * 2 + 1), WithinAbs(0.0, 1e-12));
                    REQUIRE_THAT(du_dw(1, col * 2 + 0), WithinAbs(0.0, 1e-12));
                }
            }
        }
    }

    SECTION("Differentiation consistency with interpolation")
    {
        JacobiPolynomialTpl<double, 4> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Create test data
        Eigen::Matrix<double, 3, 4> w;
        w << 1.0, 2.0, 3.0, 4.0,
             0.5, 1.5, 2.5, 3.5,
             2.0, 4.0, 6.0, 8.0;

        double t = 0.3; // Test point

        Eigen::Vector3d u;
        poly.barycentricInterpolation(t, w, u);

        Eigen::Matrix<double, 3, 12> du_dw; // 3 x (3*4)
        poly.barycentricInterpolationDiff(t, w, du_dw);

        // Check that du_dw * vec(w) ≈ u (finite difference validation)
        Eigen::VectorXd w_vec(12);
        for (int i = 0; i < 4; ++i)
        {
            w_vec.segment<3>(i * 3) = w.col(i);
        }

        Eigen::Vector3d u_reconstructed = du_dw * w_vec;

        for (int i = 0; i < 3; ++i)
        {
            REQUIRE_THAT(u_reconstructed[i], WithinAbs(u[i], 1e-10));
        }
    }

    SECTION("Differentiation matrix properties")
    {
        JacobiPolynomialTpl<double, 3> poly(0.0, 0.0);

        Eigen::Matrix<double, 2, 3> w;
        w << 1.0, 2.0, 3.0,
             4.0, 5.0, 6.0;

        double t = 0.4;
        Eigen::Matrix<double, 2, 6> du_dw;
        poly.barycentricInterpolationDiff(t, w, du_dw);

        // Each row should sum to 1 (partition of unity property)
        for (int row = 0; row < 2; ++row)
        {
            double row_sum = 0.0;
            for (int col = 0; col < 3; ++col)
            {
                row_sum += du_dw(row, col * 2 + row); // Only diagonal elements contribute
            }
            REQUIRE_THAT(row_sum, WithinAbs(1.0, 1e-10));
        }
    }
}

TEST_CASE("JacobiPolynomial - Template Parameters and Precision", "[polynomial]")
{
    SECTION("Different scalar types")
    {
        // Float precision
        JacobiPolynomialTpl<float, 3> float_poly(0.0f, 0.0f);

        // Double precision
        JacobiPolynomialTpl<double, 3> double_poly(0.0, 0.0);

        // Both should produce valid results
        const auto& float_nodes = float_poly.get_nodes();
        const auto& double_nodes = double_poly.get_nodes();

        REQUIRE(float_nodes.size() == 3);
        REQUIRE(double_nodes.size() == 3);

        // Results should be consistent within float precision
        for (int i = 0; i < 3; ++i)
        {
            REQUIRE_THAT(static_cast<double>(float_nodes[i]), WithinAbs(double_nodes[i], 1e-6));
        }
    }

    SECTION("Different polynomial degrees")
    {
        // Test various degrees
        JacobiPolynomialTpl<double, 1> poly1(0.0, 0.0);
        JacobiPolynomialTpl<double, 2> poly2(0.0, 0.0);
        JacobiPolynomialTpl<double, 5> poly5(0.0, 0.0);
        JacobiPolynomialTpl<double, 10> poly10(0.0, 0.0);

        REQUIRE(poly1.get_nodes().size() == 1);
        REQUIRE(poly2.get_nodes().size() == 2);
        REQUIRE(poly5.get_nodes().size() == 5);
        REQUIRE(poly10.get_nodes().size() == 10);
    }

    SECTION("Different Eigen storage options")
    {
        JacobiPolynomialTpl<double, 4, Eigen::ColMajor> col_major_poly(0.0, 0.0);
        JacobiPolynomialTpl<double, 4, Eigen::RowMajor> row_major_poly(0.0, 0.0);

        const auto& col_nodes = col_major_poly.get_nodes();
        const auto& row_nodes = row_major_poly.get_nodes();

        REQUIRE(col_nodes.size() == row_nodes.size());

        // Results should be identical regardless of storage order
        for (int i = 0; i < col_nodes.size(); ++i)
        {
            REQUIRE_THAT(col_nodes[i], WithinAbs(row_nodes[i], 1e-15));
        }
    }
}

TEST_CASE("JacobiPolynomial - Edge Cases and Error Conditions", "[polynomial]")
{
    SECTION("Interpolation at extreme parameter values")
    {
        // Test with large alpha, beta values
        JacobiPolynomialTpl<double, 3> extreme_poly(10.0, 15.0);
        const auto& nodes = extreme_poly.get_nodes();

        // Nodes should still be in [0,1] and ordered
        for (int i = 0; i < nodes.size(); ++i)
        {
            REQUIRE(nodes[i] >= 0.0);
            REQUIRE(nodes[i] <= 1.0);
        }

        for (int i = 0; i < nodes.size() - 1; ++i)
        {
            REQUIRE(nodes[i] < nodes[i + 1]);
        }
    }

    SECTION("Interpolation with very small values")
    {
        JacobiPolynomialTpl<double, 4> poly(0.0, 0.0);

        Eigen::Matrix<double, 1, 4> w;
        w.setConstant(1e-15);

        Eigen::Vector1d result;
        poly.barycentricInterpolation(0.5, w, result);

        REQUIRE(std::isfinite(result[0]));
        REQUIRE_THAT(result[0], WithinAbs(1e-15, 1e-16));
    }

    SECTION("Interpolation with very large values")
    {
        JacobiPolynomialTpl<double, 4> poly(0.0, 0.0);

        Eigen::Matrix<double, 1, 4> w;
        w.setConstant(1e15);

        Eigen::Vector1d result;
        poly.barycentricInterpolation(0.5, w, result);

        REQUIRE(std::isfinite(result[0]));
        REQUIRE_THAT(result[0], WithinAbs(1e15, 1e14));
    }

    SECTION("Interpolation at boundary points")
    {
        JacobiPolynomialTpl<double, 5> poly(0.0, 0.0);

        Eigen::Matrix<double, 2, 5> w;
        for (int i = 0; i < 5; ++i)
        {
            w(0, i) = i + 1.0;
            w(1, i) = (i + 1.0) * 2.0;
        }

        Eigen::Vector2d result;

        // Test at t = 0 and t = 1
        poly.barycentricInterpolation(0.0, w, result);
        REQUIRE(std::isfinite(result[0]));
        REQUIRE(std::isfinite(result[1]));

        poly.barycentricInterpolation(1.0, w, result);
        REQUIRE(std::isfinite(result[0]));
        REQUIRE(std::isfinite(result[1]));
    }

    SECTION("Zero-valued interpolation data")
    {
        JacobiPolynomialTpl<double, 3> poly(0.0, 0.0);

        Eigen::Matrix<double, 2, 3> w;
        w.setZero();

        Eigen::Vector2d result;
        poly.barycentricInterpolation(0.3, w, result);

        REQUIRE_THAT(result[0], WithinAbs(0.0, 1e-15));
        REQUIRE_THAT(result[1], WithinAbs(0.0, 1e-15));
    }
}

TEST_CASE("JacobiPolynomial - Mathematical Properties", "[polynomial]")
{
    SECTION("Partition of unity property")
    {
        JacobiPolynomialTpl<double, 5> poly(0.0, 0.0);

        // Constant function should interpolate to constant
        Eigen::Matrix<double, 1, 5> ones;
        ones.setOnes();

        Eigen::Vector1d result;

        // Test at multiple points
        for (double t = 0.0; t <= 1.0; t += 0.1)
        {
            poly.barycentricInterpolation(t, ones, result);
            REQUIRE_THAT(result[0], WithinAbs(1.0, 1e-12));
        }
    }

    SECTION("Linear reproduction property")
    {
        JacobiPolynomialTpl<double, 6> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Test multiple linear functions
        std::vector<std::pair<double, double>> linear_funcs = {
            {0.0, 1.0},   // f(x) = 1
            {1.0, 0.0},   // f(x) = x
            {2.0, 3.0},   // f(x) = 2x + 3
            {-1.0, 5.0}   // f(x) = -x + 5
        };

        for (const auto& [slope, intercept] : linear_funcs)
        {
            Eigen::Matrix<double, 1, 6> w;
            for (int i = 0; i < 6; ++i)
            {
                w(0, i) = slope * nodes[i] + intercept;
            }

            Eigen::Vector1d result;

            // Test at multiple points
            for (double t = 0.1; t <= 0.9; t += 0.2)
            {
                poly.barycentricInterpolation(t, w, result);
                double expected = slope * t + intercept;
                REQUIRE_THAT(result[0], WithinAbs(expected, 1e-11));
            }
        }
    }

    SECTION("Symmetry properties for Legendre polynomials")
    {
        JacobiPolynomialTpl<double, 6> legendre(0.0, 0.0); // Legendre case
        const auto& nodes = legendre.get_nodes();

        // For Legendre polynomials, nodes should be symmetric about 0.5
        for (int i = 0; i < nodes.size(); ++i)
        {
            double symmetric_node = 1.0 - nodes[nodes.size() - 1 - i];
            REQUIRE_THAT(nodes[i], WithinAbs(symmetric_node, 1e-12));
        }
    }

    SECTION("Interpolation exactness for polynomial degree N-1")
    {
        JacobiPolynomialTpl<double, 5> poly(0.0, 0.0);
        const auto& nodes = poly.get_nodes();

        // Test with a polynomial of degree 4 (N-1)
        // p(x) = x^4 - 2*x^3 + 3*x^2 - x + 2
        auto polynomial = [](double x) {
            return x*x*x*x - 2.0*x*x*x + 3.0*x*x - x + 2.0;
        };

        Eigen::Matrix<double, 1, 5> w;
        for (int i = 0; i < 5; ++i)
        {
            w(0, i) = polynomial(nodes[i]);
        }

        Eigen::Vector1d result;

        // Should interpolate exactly at any point
        for (double t = 0.05; t <= 0.95; t += 0.15)
        {
            poly.barycentricInterpolation(t, w, result);
            double expected = polynomial(t);
            REQUIRE_THAT(result[0], WithinAbs(expected, 1e-10));
        }
    }
}
