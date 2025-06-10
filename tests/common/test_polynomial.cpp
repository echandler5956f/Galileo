#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include "galileo/common/polynomial/polynomial.hpp"
#include "../helpers/catch_eigen_matchers.hpp"

using namespace galileo;
using namespace GalileoMatchers;

TEMPLATE_TEST_CASE("Jacobi Polynomial Interpolation", "[common][polynomial]", double)
{
    using Scalar = TestType;
    const int NOrder = 4; // Using a 4th order polynomial for a non-trivial test
    const int UDim = 2;   // Dimension of the interpolated vector

    using Poly_t = math::JacobiPolynomialTpl<Scalar, NOrder, 0>;
    using MatrixU_N = Eigen::Matrix<Scalar, UDim, NOrder>;
    using VectorU = Eigen::Matrix<Scalar, UDim, 1>;
    using MatrixU_UN = Eigen::Matrix<Scalar, UDim, UDim * NOrder>;

    Poly_t poly(0., 0.); // Legendre polynomial
    auto nodes = poly.get_nodes();

    // A known quadratic function to test against: u(t) = [t^2, 2t]
    auto test_func = [](Scalar t) -> VectorU
    {
        return VectorU(t * t, 2 * t);
    };

    // Evaluate the function at the nodes to get the control points
    MatrixU_N w;
    for (int i = 0; i < NOrder; ++i)
    {
        w.col(i) = test_func(nodes(i));
    }

    SECTION("Barycentric Interpolation Evaluation")
    {
        Scalar t = 0.3; // A point to interpolate

        VectorU u_interpolated;
        poly.barycentricInterpolation(t, w, u_interpolated);

        VectorU u_analytical = test_func(t);
        REQUIRE_THAT(u_interpolated, Approx(u_analytical, 1e-7));
    }

    SECTION("Barycentric Interpolation Derivative Verification")
    {
        Scalar t = 0.6; // A different point to test derivative

        MatrixU_UN du_dw_analytical;
        poly.barycentricInterpolationDiff(t, w, du_dw_analytical);

        MatrixU_UN du_dw_numerical = MatrixU_UN::Zero();
        Scalar epsilon = 1e-7;

        for (int i = 0; i < w.size(); ++i)
        {
            MatrixU_N w_plus = w;
            w_plus(i) += epsilon;
            VectorU u_plus;
            poly.barycentricInterpolation(t, w_plus, u_plus);

            MatrixU_N w_minus = w;
            w_minus(i) -= epsilon;
            VectorU u_minus;
            poly.barycentricInterpolation(t, w_minus, u_minus);

            // The i-th column of the numerical derivative matrix
            // Note: Eigen matrices are column-major by default. w(i) accesses elements
            // column by column. The derivative df/dw_ij is placed in column (j*rows + i)
            // for a matrix w of size (rows x cols). This logic is correct for mapping
            // a matrix perturbation to a column in the Jacobian.
            du_dw_numerical.col(i) = (u_plus - u_minus) / (2 * epsilon);
        }
        REQUIRE_THAT(du_dw_analytical, Approx(du_dw_numerical, epsilon * 100));
    }
} 