#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/core/states/implementations/euclidean.hpp"
#include "galileo/core/basic-spec.hpp"
#include "galileo/multibody/robot-spec.hpp"
#include "galileo/core/states/state-base.hpp"
#include "galileo/common/meta/dimension.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

using namespace galileo;
using namespace Catch::Matchers;

template <typename Scalar>
constexpr Scalar TOLERANCE = std::numeric_limits<Scalar>::epsilon() * 1000;

namespace galileo {

// Forward declarations for dummy actuation classes
template <typename RobotSpec>
class DummyActuationTpl;

template <typename RobotSpec>
struct DummyActuationDataTpl;

// Dummy actuation classes for testing
template <typename RobotSpec>
class DummyActuationTpl
{
public:
    using RS = RobotSpec;
    using Data_t = DummyActuationDataTpl<RS>;
};

template <typename RobotSpec>
struct DummyActuationDataTpl
{
    using RS = RobotSpec;
};

// Specialization of traits for dummy actuation
template <typename RobotSpec>
struct traits<DummyActuationTpl<RobotSpec>>
{
    using Model_t = DummyActuationTpl<RobotSpec>;
    using Data_t = DummyActuationDataTpl<RobotSpec>;
};

} // namespace galileo

// Test helper functions
template <typename StateType>
void test_basic_state_properties(StateType& state, int expected_nx, int expected_ndx, int expected_nq, int expected_nv)
{
    REQUIRE(state.get_nx() == expected_nx);
    REQUIRE(state.get_ndx() == expected_ndx);
    REQUIRE(state.get_nq() == expected_nq);
    REQUIRE(state.get_nv() == expected_nv);

    // Test zero state
    typename StateType::VectorNx_t zero_state = state.zero();
    REQUIRE(zero_state.size() == expected_nx);
    for (int i = 0; i < expected_nx; ++i) {
        REQUIRE(zero_state(i) == 0.0);
    }

    // Test random state
    typename StateType::VectorNx_t rand_state = state.rand();
    REQUIRE(rand_state.size() == expected_nx);
}

template <typename StateType>
void test_manifold_operations(StateType& state, int nx, int ndx)
{
    // Create test vectors
    typename StateType::VectorNx_t x0 = StateType::VectorNx_t::Random(nx);
    typename StateType::VectorNx_t x1 = StateType::VectorNx_t::Random(nx);
    typename StateType::VectorNdx_t dx = StateType::VectorNdx_t::Random(ndx);

    // Test diff
    typename StateType::VectorNdx_t dx_out = StateType::VectorNdx_t::Zero(ndx);

    state.diff(x0, x1, dx_out);

    REQUIRE(dx_out.size() == ndx);
    for (int i = 0; i < ndx; ++i) {
        REQUIRE_THAT(dx_out(i), WithinAbs(x1(i) - x0(i), TOLERANCE<double>));
    }

    // Test integrate
    typename StateType::VectorNx_t x_out(nx);
    state.integrate(x0, dx, x_out);

    REQUIRE(x_out.size() == nx);
    for (int i = 0; i < ndx; ++i) {
        REQUIRE_THAT(x_out(i), WithinAbs(x0(i) + dx(i), TOLERANCE<double>));
    }

    // Test consistency: integrate(x0, diff(x0, x1)) == x1
    typename StateType::VectorNdx_t dx_test = state.diff_dx(x0, x1);
    typename StateType::VectorNx_t x1_reconstructed = state.integrate_x(x0, dx_test);

    REQUIRE(x1_reconstructed.size() == nx);
    for (int i = 0; i < nx; ++i) {
        REQUIRE_THAT(x1_reconstructed(i), WithinAbs(x1(i), TOLERANCE<double>));
    }
}

TEST_CASE("StateEuclideanTpl - Fixed Dimensions", "[euclidean][fixed]")
{
    SECTION("Standard robot configuration")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 6, 6, 6, 6, 0, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(24, -10.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(24, 10.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test dimension accessors
        REQUIRE(state.get_nqb() == 6);
        REQUIRE(state.get_nqj() == 6);
        REQUIRE(state.get_nq() == 12);
        REQUIRE(state.get_nvb() == 6);
        REQUIRE(state.get_nvj() == 6);
        REQUIRE(state.get_nv() == 12);
        REQUIRE(state.get_nrotors() == 0);
        REQUIRE(state.get_nx() == 24);
        REQUIRE(state.get_ndx() == 24);
        REQUIRE(state.get_nua() == 6);

        // Test bounds
        REQUIRE(state.get_lb().size() == 24);
        REQUIRE(state.get_ub().size() == 24);
        for (int i = 0; i < 24; ++i) {
            REQUIRE(state.get_lb()(i) == -10.0);
            REQUIRE(state.get_ub()(i) == 10.0);
        }

        // Test basic state properties
        test_basic_state_properties(state, 24, 24, 12, 12);

        // Test manifold operations
        test_manifold_operations(state, 24, 24);
    }

    SECTION("No floating base")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 0, 6, 0, 6, 0, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(12, -5.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(12, 5.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test dimension accessors
        REQUIRE(state.get_nqb() == 0);
        REQUIRE(state.get_nqj() == 6);
        REQUIRE(state.get_nq() == 6);
        REQUIRE(state.get_nvb() == 0);
        REQUIRE(state.get_nvj() == 6);
        REQUIRE(state.get_nv() == 6);
        REQUIRE(state.get_nrotors() == 0);
        REQUIRE(state.get_nx() == 12);
        REQUIRE(state.get_ndx() == 12);
        REQUIRE(state.get_nua() == 6);

        // Test basic state properties
        test_basic_state_properties(state, 12, 12, 6, 6);

        // Test manifold operations
        test_manifold_operations(state, 12, 12);
    }

    SECTION("Minimal configuration")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 1, 1, 1, 1, 1, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(4, -1.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(4, 1.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test dimension accessors
        REQUIRE(state.get_nqb() == 1);
        REQUIRE(state.get_nqj() == 1);
        REQUIRE(state.get_nq() == 2);
        REQUIRE(state.get_nvb() == 1);
        REQUIRE(state.get_nvj() == 1);
        REQUIRE(state.get_nv() == 2);
        REQUIRE(state.get_nrotors() == 1);
        REQUIRE(state.get_nx() == 4);
        REQUIRE(state.get_ndx() == 4);
        REQUIRE(state.get_nua() == 2);

        // Test basic state properties
        test_basic_state_properties(state, 4, 4, 2, 2);

        // Test manifold operations
        test_manifold_operations(state, 4, 4);
    }
}

TEST_CASE("StateEuclideanTpl - Dynamic Dimensions", "[euclidean][dynamic]")
{
    SECTION("Dynamic robot configuration")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;

        // Setup dynamic dimensions
        rs.NQb_dim.set_value(6);
        rs.NQj_dim.set_value(6);
        rs.NVb_dim.set_value(6);
        rs.NVj_dim.set_value(6);
        rs.NRotors_dim.set_value(0);
        rs.NQ_dim.set_value(12);
        rs.NV_dim.set_value(12);
        rs.NX_dim.set_value(24);
        rs.NDX_dim.set_value(24);
        rs.NUa_dim.set_value(6);

        RS::VectorNx_t lb(24);
        RS::VectorNx_t ub(24);
        lb.setConstant(-8.0);
        ub.setConstant(8.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test dimension accessors
        REQUIRE(state.get_nqb() == 6);
        REQUIRE(state.get_nqj() == 6);
        REQUIRE(state.get_nq() == 12);
        REQUIRE(state.get_nvb() == 6);
        REQUIRE(state.get_nvj() == 6);
        REQUIRE(state.get_nv() == 12);
        REQUIRE(state.get_nrotors() == 0);
        REQUIRE(state.get_nx() == 24);
        REQUIRE(state.get_ndx() == 24);
        REQUIRE(state.get_nua() == 6);

        // Test bounds
        REQUIRE(state.get_lb().size() == 24);
        REQUIRE(state.get_ub().size() == 24);
        for (int i = 0; i < 24; ++i) {
            REQUIRE(state.get_lb()(i) == -8.0);
            REQUIRE(state.get_ub()(i) == 8.0);
        }

        // Test basic state properties
        test_basic_state_properties(state, 24, 24, 12, 12);

        // Test manifold operations
        test_manifold_operations(state, 24, 24);
    }
}

TEST_CASE("StateEuclideanTpl - Jacobian Operations", "[euclidean][jacobians]")
{
    SECTION("Jdiff - Fixed dimensions")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 6, 6, 6, 6, 0, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(24, -10.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(24, 10.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Create test vectors
        typename RS::VectorNx_t x0 = RS::VectorNx_t::Random(24);
        typename RS::VectorNx_t x1 = RS::VectorNx_t::Random(24);

        // Test Jdiff with both components
        RS::MatrixNdx_t Jfirst = RS::MatrixNdx_t::Zero(24, 24);
        RS::MatrixNdx_t Jsecond = RS::MatrixNdx_t::Zero(24, 24);

        state.Jdiff(x0, x1, Jfirst, Jsecond, both);

        // For Euclidean space, Jdiff should be identity matrices with appropriate signs
        for (int i = 0; i < 24; ++i) {
            for (int j = 0; j < 24; ++j) {
                if (i == j) {
                    REQUIRE_THAT(Jfirst(i, j), WithinAbs(-1.0, TOLERANCE<double>));
                    REQUIRE_THAT(Jsecond(i, j), WithinAbs(1.0, TOLERANCE<double>));
                } else {
                    REQUIRE_THAT(Jfirst(i, j), WithinAbs(0.0, TOLERANCE<double>));
                    REQUIRE_THAT(Jsecond(i, j), WithinAbs(0.0, TOLERANCE<double>));
                }
            }
        }

        // Test Jdiff with first component only
        RS::MatrixNdx_t Jfirst_only = RS::MatrixNdx_t::Zero(24, 24);
        RS::MatrixNdx_t Jsecond_only = RS::MatrixNdx_t::Zero(24, 24);

        state.Jdiff(x0, x1, Jfirst_only, Jsecond_only, first);

        for (int i = 0; i < 24; ++i) {
            for (int j = 0; j < 24; ++j) {
                if (i == j) {
                    REQUIRE_THAT(Jfirst_only(i, j), WithinAbs(-1.0, TOLERANCE<double>));
                } else {
                    REQUIRE_THAT(Jfirst_only(i, j), WithinAbs(0.0, TOLERANCE<double>));
                }
                // Jsecond_only should remain zero
                REQUIRE_THAT(Jsecond_only(i, j), WithinAbs(0.0, TOLERANCE<double>));
            }
        }

        // Test Jdiff_Js convenience method
        std::vector<typename RS::MatrixNdx_t> Jacs = state.Jdiff_Js(x0, x1, both);
        REQUIRE(Jacs.size() == 2);
        REQUIRE(Jacs[0].rows() == 24);
        REQUIRE(Jacs[0].cols() == 24);
        REQUIRE(Jacs[1].rows() == 24);
        REQUIRE(Jacs[1].cols() == 24);
    }

    SECTION("Jintegrate - Fixed dimensions")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 1, 1, 1, 1, 1, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(4, -10.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(4, 10.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Create test vectors
        typename RS::VectorNx_t x = RS::VectorNx_t::Random(4);
        typename RS::VectorNdx_t dx = RS::VectorNdx_t::Random(4);

        // Test Jintegrate with setto operation
        RS::MatrixNdx_t Jfirst = RS::MatrixNdx_t::Zero(4, 4);
        RS::MatrixNdx_t Jsecond = RS::MatrixNdx_t::Zero(4, 4);

        state.Jintegrate(x, dx, Jfirst, Jsecond, both, setto);

        // For Euclidean space, Jintegrate should be identity matrices
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i == j) {
                    REQUIRE_THAT(Jfirst(i, j), WithinAbs(1.0, TOLERANCE<double>));
                    REQUIRE_THAT(Jsecond(i, j), WithinAbs(1.0, TOLERANCE<double>));
                } else {
                    REQUIRE_THAT(Jfirst(i, j), WithinAbs(0.0, TOLERANCE<double>));
                    REQUIRE_THAT(Jsecond(i, j), WithinAbs(0.0, TOLERANCE<double>));
                }
            }
        }

        // Test Jintegrate with addto operation
        RS::MatrixNdx_t Jfirst_add = RS::MatrixNdx_t::Ones(4, 4);
        RS::MatrixNdx_t Jsecond_add = RS::MatrixNdx_t::Ones(4, 4);

        state.Jintegrate(x, dx, Jfirst_add, Jsecond_add, both, addto);

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i == j) {
                    REQUIRE_THAT(Jfirst_add(i, j), WithinAbs(2.0, TOLERANCE<double>));
                    REQUIRE_THAT(Jsecond_add(i, j), WithinAbs(2.0, TOLERANCE<double>));
                } else {
                    REQUIRE_THAT(Jfirst_add(i, j), WithinAbs(1.0, TOLERANCE<double>));
                    REQUIRE_THAT(Jsecond_add(i, j), WithinAbs(1.0, TOLERANCE<double>));
                }
            }
        }

        // Test Jintegrate with rmfrom operation
        RS::MatrixNdx_t Jfirst_sub = RS::MatrixNdx_t::Ones(4, 4);
        RS::MatrixNdx_t Jsecond_sub = RS::MatrixNdx_t::Ones(4, 4);

        state.Jintegrate(x, dx, Jfirst_sub, Jsecond_sub, both, rmfrom);

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i == j) {
                    REQUIRE_THAT(Jfirst_sub(i, j), WithinAbs(0.0, TOLERANCE<double>));
                    REQUIRE_THAT(Jsecond_sub(i, j), WithinAbs(0.0, TOLERANCE<double>));
                } else {
                    REQUIRE_THAT(Jfirst_sub(i, j), WithinAbs(1.0, TOLERANCE<double>));
                    REQUIRE_THAT(Jsecond_sub(i, j), WithinAbs(1.0, TOLERANCE<double>));
                }
            }
        }

        // Test Jintegrate_Js convenience method
        std::vector<typename RS::MatrixNdx_t> Jacs = state.Jintegrate_Js(x, dx, both);
        REQUIRE(Jacs.size() == 2);
        REQUIRE(Jacs[0].rows() == 4);
        REQUIRE(Jacs[0].cols() == 4);
        REQUIRE(Jacs[1].rows() == 4);
        REQUIRE(Jacs[1].cols() == 4);
    }

    SECTION("JintegrateTransport - Euclidean is trivial")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 1, 1, 1, 1, 1, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(4, -10.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(4, 10.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Create test vectors
        typename RS::VectorNx_t x = RS::VectorNx_t::Random(4);
        typename RS::VectorNdx_t dx = RS::VectorNdx_t::Random(4);
        RS::MatrixNdx_t Jin = RS::MatrixNdx_t::Random(4, 4);
        RS::MatrixNdx_t Jin_copy = Jin;

        // For Euclidean space, JintegrateTransport should do nothing
        state.JintegrateTransport(x, dx, Jin, first);

        // Matrix should be unchanged
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                REQUIRE_THAT(Jin(i, j), WithinAbs(Jin_copy(i, j), TOLERANCE<double>));
            }
        }
    }

    SECTION("Dynamic dimensions Jacobians")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;

        // Setup dynamic dimensions
        rs.NQb_dim.set_value(2);
        rs.NQj_dim.set_value(2);
        rs.NVb_dim.set_value(2);
        rs.NVj_dim.set_value(2);
        rs.NRotors_dim.set_value(1);
        rs.NQ_dim.set_value(4);
        rs.NV_dim.set_value(4);
        rs.NX_dim.set_value(8);
        rs.NDX_dim.set_value(8);
        rs.NUa_dim.set_value(3);

        RS::VectorNx_t lb(8);
        RS::VectorNx_t ub(8);
        lb.setConstant(-10.0);
        ub.setConstant(10.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Create test vectors
        typename RS::VectorNx_t x0 = RS::VectorNx_t::Random(8);
        typename RS::VectorNx_t x1 = RS::VectorNx_t::Random(8);
        typename RS::VectorNx_t x = RS::VectorNx_t::Random(8);
        typename RS::VectorNdx_t dx = RS::VectorNdx_t::Random(8);

        // Test Jdiff
        RS::MatrixNdx_t Jfirst = RS::MatrixNdx_t::Zero(8, 8);
        RS::MatrixNdx_t Jsecond = RS::MatrixNdx_t::Zero(8, 8);

        state.Jdiff(x0, x1, Jfirst, Jsecond, both);

        // Should still be identity matrices with appropriate signs
        for (int i = 0; i < 8; ++i) {
            REQUIRE_THAT(Jfirst(i, i), WithinAbs(-1.0, TOLERANCE<double>));
            REQUIRE_THAT(Jsecond(i, i), WithinAbs(1.0, TOLERANCE<double>));
        }

        // Test Jintegrate
        RS::MatrixNdx_t Jfirst_int = RS::MatrixNdx_t::Zero(8, 8);
        RS::MatrixNdx_t Jsecond_int = RS::MatrixNdx_t::Zero(8, 8);

        state.Jintegrate(x, dx, Jfirst_int, Jsecond_int, both, setto);

        // Should be identity matrices
        for (int i = 0; i < 8; ++i) {
            REQUIRE_THAT(Jfirst_int(i, i), WithinAbs(1.0, TOLERANCE<double>));
            REQUIRE_THAT(Jsecond_int(i, i), WithinAbs(1.0, TOLERANCE<double>));
        }
    }
}

TEST_CASE("StateEuclideanTpl - Bounds Management", "[euclidean][bounds]")
{
    SECTION("Bounds modification")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 6, 6, 6, 6, 0, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(24, -10.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(24, 10.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test initial bounds
        for (int i = 0; i < 24; ++i) {
            REQUIRE(state.get_lb()(i) == -10.0);
            REQUIRE(state.get_ub()(i) == 10.0);
        }

        // Test bounds modification
        RS::VectorNx_t new_lb = RS::VectorNx_t::Constant(24, -5.0);
        RS::VectorNx_t new_ub = RS::VectorNx_t::Constant(24, 15.0);

        state.set_lb(new_lb);
        state.set_ub(new_ub);

        // Test modified bounds
        for (int i = 0; i < 24; ++i) {
            REQUIRE(state.get_lb()(i) == -5.0);
            REQUIRE(state.get_ub()(i) == 15.0);
        }
    }

    SECTION("Bounds with different values")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 1, 1, 1, 1, 1, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb(4);
        RS::VectorNx_t ub(4);

        lb << -1.0, -2.0, -3.0, -4.0;
        ub << 1.0, 2.0, 3.0, 4.0;

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test heterogeneous bounds
        REQUIRE(state.get_lb()(0) == -1.0);
        REQUIRE(state.get_lb()(1) == -2.0);
        REQUIRE(state.get_lb()(2) == -3.0);
        REQUIRE(state.get_lb()(3) == -4.0);

        REQUIRE(state.get_ub()(0) == 1.0);
        REQUIRE(state.get_ub()(1) == 2.0);
        REQUIRE(state.get_ub()(2) == 3.0);
        REQUIRE(state.get_ub()(3) == 4.0);
    }

    SECTION("Dynamic dimensions bounds")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;

        // Setup dynamic dimensions
        rs.NQb_dim.set_value(1);
        rs.NQj_dim.set_value(1);
        rs.NVb_dim.set_value(1);
        rs.NVj_dim.set_value(1);
        rs.NRotors_dim.set_value(1);
        rs.NQ_dim.set_value(2);
        rs.NV_dim.set_value(2);
        rs.NX_dim.set_value(4);
        rs.NDX_dim.set_value(4);
        rs.NUa_dim.set_value(2);

        RS::VectorNx_t lb(4);
        RS::VectorNx_t ub(4);
        lb.setConstant(-20.0);
        ub.setConstant(20.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test bounds
        for (int i = 0; i < 4; ++i) {
            REQUIRE(state.get_lb()(i) == -20.0);
            REQUIRE(state.get_ub()(i) == 20.0);
        }

        // Test bounds modification
        RS::VectorNx_t new_lb(4);
        RS::VectorNx_t new_ub(4);

        new_lb << -10.0, -15.0, -20.0, -25.0;
        new_ub << 10.0, 15.0, 20.0, 25.0;

        state.set_lb(new_lb);
        state.set_ub(new_ub);

        REQUIRE(state.get_lb()(0) == -10.0);
        REQUIRE(state.get_lb()(1) == -15.0);
        REQUIRE(state.get_lb()(2) == -20.0);
        REQUIRE(state.get_lb()(3) == -25.0);

        REQUIRE(state.get_ub()(0) == 10.0);
        REQUIRE(state.get_ub()(1) == 15.0);
        REQUIRE(state.get_ub()(2) == 20.0);
        REQUIRE(state.get_ub()(3) == 25.0);
    }
}

TEST_CASE("StateEuclideanTpl - Edge Cases", "[euclidean][edge_cases]")
{
    SECTION("Zero state dimension")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 0, 0, 0, 0, 0, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(0, -10.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(0, 10.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test dimension accessors
        REQUIRE(state.get_nqb() == 0);
        REQUIRE(state.get_nqj() == 0);
        REQUIRE(state.get_nq() == 0);
        REQUIRE(state.get_nvb() == 0);
        REQUIRE(state.get_nvj() == 0);
        REQUIRE(state.get_nv() == 0);
        REQUIRE(state.get_nrotors() == 0);
        REQUIRE(state.get_nx() == 0);
        REQUIRE(state.get_ndx() == 0);
        REQUIRE(state.get_nua() == 0);

        // Test zero state
        typename RS::VectorNx_t zero_state = state.zero();
        REQUIRE(zero_state.size() == 0);

        // Test random state
        typename RS::VectorNx_t rand_state = state.rand();
        REQUIRE(rand_state.size() == 0);
    }

    SECTION("Large state dimension")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 50, 50, 50, 50, 10, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(200, -100.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(200, 100.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test dimension accessors
        REQUIRE(state.get_nqb() == 50);
        REQUIRE(state.get_nqj() == 50);
        REQUIRE(state.get_nq() == 100);
        REQUIRE(state.get_nvb() == 50);
        REQUIRE(state.get_nvj() == 50);
        REQUIRE(state.get_nv() == 100);
        REQUIRE(state.get_nrotors() == 10);
        REQUIRE(state.get_nx() == 200);
        REQUIRE(state.get_ndx() == 200);
        REQUIRE(state.get_nua() == 60);

        // Test basic state properties
        test_basic_state_properties(state, 200, 200, 100, 100);
    }

    SECTION("Random state generation consistency")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 2, 2, 2, 2, 1, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        RS::VectorNx_t lb = RS::VectorNx_t::Constant(8, -1.0);
        RS::VectorNx_t ub = RS::VectorNx_t::Constant(8, 1.0);

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Generate multiple random states and check they're different
        typename RS::VectorNx_t rand_state1 = state.rand();
        typename RS::VectorNx_t rand_state2 = state.rand();
        typename RS::VectorNx_t rand_state3 = state.rand();

        REQUIRE(rand_state1.size() == 8);
        REQUIRE(rand_state2.size() == 8);
        REQUIRE(rand_state3.size() == 8);

        // Check that at least some values are different
        bool different_12 = false, different_13 = false, different_23 = false;
        for (int i = 0; i < 8; ++i) {
            if (std::abs(rand_state1(i) - rand_state2(i)) > TOLERANCE<double>) {
                different_12 = true;
            }
            if (std::abs(rand_state1(i) - rand_state3(i)) > TOLERANCE<double>) {
                different_13 = true;
            }
            if (std::abs(rand_state2(i) - rand_state3(i)) > TOLERANCE<double>) {
                different_23 = true;
            }
        }

        REQUIRE(different_12);
        REQUIRE(different_13);
        REQUIRE(different_23);
    }
}

TEMPLATE_TEST_CASE("StateEuclideanTpl - Different Scalar Types", "[euclidean][scalar_types]", float, double, long double)
{
    using Scalar = TestType;
    constexpr Scalar tolerance = TOLERANCE<Scalar>;

    SECTION("Basic operations with different scalar types")
    {
        using BS = BasicSpecTpl<Scalar, Scalar, 0>;
        using RS = RobotSpecTpl<BS, 2, 2, 2, 2, 1, StateEuclideanTpl, DummyActuationTpl>;

        RS rs;
        typename RS::VectorNx_t lb = RS::VectorNx_t::Constant(8, Scalar(-5.0));
        typename RS::VectorNx_t ub = RS::VectorNx_t::Constant(8, Scalar(5.0));

        StateEuclideanTpl<RS> state(rs, lb, ub);

        // Test dimension accessors
        REQUIRE(state.get_nx() == 8);
        REQUIRE(state.get_ndx() == 8);
        REQUIRE(state.get_nq() == 4);
        REQUIRE(state.get_nv() == 4);

        // Test zero state
        typename RS::VectorNx_t zero_state = state.zero();
        REQUIRE(zero_state.size() == 8);
        for (int i = 0; i < 8; ++i) {
            REQUIRE(zero_state(i) == Scalar(0.0));
        }

        // Test manifold operations
        typename RS::VectorNx_t x0 = RS::VectorNx_t::Random(8);
        typename RS::VectorNx_t x1 = RS::VectorNx_t::Random(8);
        typename RS::VectorNdx_t dx = RS::VectorNdx_t::Random(8);

        // Test diff
        typename RS::VectorNdx_t dx_out;
        state.diff(x0, x1, dx_out);

        REQUIRE(dx_out.size() == 8);
        for (int i = 0; i < 8; ++i) {
            REQUIRE_THAT(static_cast<double>(dx_out(i)), WithinAbs(static_cast<double>(x1(i) - x0(i)), tolerance));
        }

        // Test integrate
        typename RS::VectorNx_t x_out;
        state.integrate(x0, dx, x_out);

        REQUIRE(x_out.size() == 8);
        for (int i = 0; i < 8; ++i) {
            REQUIRE_THAT(static_cast<double>(x_out(i)), WithinAbs(static_cast<double>(x0(i) + dx(i)), tolerance));
        }

        // Test consistency
        typename RS::VectorNdx_t dx_test = state.diff_dx(x0, x1);
        typename RS::VectorNx_t x1_reconstructed = state.integrate_x(x0, dx_test);

        REQUIRE(x1_reconstructed.size() == 8);
        for (int i = 0; i < 8; ++i) {
            REQUIRE_THAT(static_cast<double>(x1_reconstructed(i)), WithinAbs(static_cast<double>(x1(i)), tolerance));
        }
    }
}
