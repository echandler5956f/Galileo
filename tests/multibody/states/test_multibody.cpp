#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include "galileo/multibody/actuations/implementations/actuation-floating-base.hpp"
#include "galileo/multibody/states/implementations/state-multibody.hpp"

#include <filesystem>
#include <iostream>
#include <random>

#include "utils/resource_finder.hpp"

using namespace galileo;
using namespace Catch::Matchers;

namespace test_helpers
{

    // Robot model specifications
    struct RobotModelSpec
    {
        std::string name;
        int nqj;
        int nvj;
        bool has_floating_base;

        std::string urdf_path() const
        {
            return galileo::testing::get_robot_urdf_path(name);
        }
    };

    static const std::vector<RobotModelSpec> robot_specs = {
        {"atlas", 30, 30, true},
        {"huron", 12, 12, true},
        {"go1", 12, 12, true}};

    // Tolerance helper
    template <typename Scalar>
    constexpr Scalar tolerance()
    {
        if constexpr (std::is_same_v<Scalar, float>)
        {
            return 1e-5f;
        }
        else if constexpr (std::is_same_v<Scalar, double>)
        {
            return 1e-10;
        }
        else
        {
            return 1e-12L;
        }
    }

    // Test fixture for robot models
    template <typename RobotSpec>
    class MultibodyStateTestFixture
    {
    public:
        using RS = RobotSpec;
        using StateType = StateMultibodyTpl<RS>;
        using Scalar = typename RS::VarScalar;

        std::unique_ptr<typename RS::RobotModel_t> model;
        std::unique_ptr<typename RS::RobotData_t> data;
        std::unique_ptr<RS> rs;
        std::unique_ptr<StateType> state;

        void setup_robot(const RobotModelSpec &spec)
        {
            // Load URDF
            pinocchio::ModelTpl<double> model_tmp;
            if (spec.has_floating_base)
            {
                pinocchio::urdf::buildModel(spec.urdf_path(), pinocchio::JointModelFreeFlyer(), model_tmp);
            }
            else
            {
                pinocchio::urdf::buildModel(spec.urdf_path(), model_tmp);
            }

            pinocchio::DataTpl<double> data_tmp = pinocchio::DataTpl<double>(model_tmp);

            model = std::make_unique<typename RS::RobotModel_t>(model_tmp);
            data = std::make_unique<typename RS::RobotData_t>(data_tmp);

            rs = std::make_unique<RS>();

            // Create state
            state = std::make_unique<StateType>(*model);

            // Validate dimensions
            REQUIRE(IsValidRobotSpec(*rs));

            if (spec.has_floating_base)
            {
                REQUIRE(state->get_nqb() == 7);
                REQUIRE(state->get_nvb() == 6);
            }
            REQUIRE(state->get_nqj() == spec.nqj);
            REQUIRE(state->get_nvj() == spec.nvj);
        }

        // Generate random valid configuration
        typename RS::VectorNx_t random_state()
        {
            return state->rand();
        }

        // Generate random velocity
        typename RS::VectorNdx_t random_velocity()
        {
            return RS::VectorNdx_t::Random(state->get_ndx());
        }
    };

} // namespace test_helpers

using namespace test_helpers;

// Test dimension management for fixed and dynamic dimensions
TEST_CASE("StateMultibodyTpl - Dimension Management", "[multibody][dimensions]")
{
    SECTION("Fixed dimensions with floating base")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 7, 30, 6, 30, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

        MultibodyStateTestFixture<RS> fixture;

        RobotModelSpec spec = {"atlas", 30, 30, true};
        fixture.setup_robot(spec);

        REQUIRE(fixture.state->get_nq() == 37);
        REQUIRE(fixture.state->get_nv() == 36);
        REQUIRE(fixture.state->get_nx() == 73);
        REQUIRE(fixture.state->get_ndx() == 72);
        REQUIRE(fixture.state->get_nua() == 30);
    }

    SECTION("Dynamic dimensions")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

        for (const auto &spec : robot_specs)
        {
            MultibodyStateTestFixture<RS> fixture;
            fixture.setup_robot(spec);

            // Verify dimensions were set correctly from model
            if (spec.has_floating_base)
            {
                REQUIRE(fixture.state->get_nqb() == 7);
                REQUIRE(fixture.state->get_nvb() == 6);
                REQUIRE(fixture.state->get_nq() == 7 + spec.nqj);
                REQUIRE(fixture.state->get_nv() == 6 + spec.nvj);
            }
        }
    }

    SECTION("Mixed fixed/dynamic dimensions")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS = RobotSpecTpl<BS, 7, detail::Dynamic, 6, detail::Dynamic, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

        MultibodyStateTestFixture<RS> fixture;
        fixture.setup_robot({"go1", 12, 12, true});

        REQUIRE(fixture.state->get_nqb() == 7);
        REQUIRE(fixture.state->get_nvb() == 6);
        REQUIRE(fixture.state->get_nqj() == 12);
        REQUIRE(fixture.state->get_nvj() == 12);
    }
}

// Test state generation methods
TEST_CASE("StateMultibodyTpl - State Generation", "[multibody][state_generation]")
{
    using BS = BasicSpecTpl<double, double, 0>;
    using RS = RobotSpecTpl<BS, 7, 12, 6, 12, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

    MultibodyStateTestFixture<RS> fixture;
    fixture.setup_robot({"go1", 12, 12, true});

    SECTION("Zero state")
    {
        Eigen::GMatrix<double, 37, 1> x0 = fixture.state->zero();
        REQUIRE(x0.size() == fixture.state->get_nx());

        // Configuration should be neutral
        Eigen::VectorXd q_neutral = pinocchio::neutral(*fixture.model);
        for (int i = 0; i < fixture.state->get_nq(); ++i)
        {
            REQUIRE_THAT(x0(i), WithinAbs(q_neutral(i), tolerance<double>()));
        }

        // Velocities should be zero
        for (int i = fixture.state->get_nq(); i < fixture.state->get_nx(); ++i)
        {
            REQUIRE_THAT(x0(i), WithinAbs(0.0, tolerance<double>()));
        }
    }

    SECTION("Random state")
    {
        // Generate multiple random states and verify they're different
        Eigen::GMatrix<double, 37, 1> x1 = fixture.state->rand();
        Eigen::GMatrix<double, 37, 1> x2 = fixture.state->rand();

        REQUIRE(x1.size() == fixture.state->get_nx());
        REQUIRE(x2.size() == fixture.state->get_nx());

        // Check states are different
        bool different = false;
        for (int i = 0; i < fixture.state->get_nx(); ++i)
        {
            if (std::abs(x1(i) - x2(i)) > tolerance<double>())
            {
                different = true;
                break;
            }
        }
        REQUIRE(different);

        // Configuration should be valid
        Eigen::GMatrix<double, 19, 1> q1 = x1.head(fixture.state->get_nq());
        REQUIRE(pinocchio::isNormalized(*fixture.model, q1, tolerance<double>() * 100));
    }
}

// Test manifold operations
TEST_CASE("StateMultibodyTpl - Manifold Operations", "[multibody][manifold]")
{
    using BS = BasicSpecTpl<double, double, 0>;
    using RS = RobotSpecTpl<BS, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

    for (const auto &spec : robot_specs)
    {
        SECTION(spec.name + " manifold operations")
        {
            MultibodyStateTestFixture<RS> fixture;
            fixture.setup_robot(spec);

            typename RS::VectorNx_t x0 = fixture.random_state();
            typename RS::VectorNx_t x1 = fixture.random_state();

            SECTION("Diff operation")
            {
                typename RS::VectorNdx_t dx(fixture.state->get_ndx());
                fixture.state->diff(x0, x1, dx);

                // Verify configuration difference
                typename RS::VectorNv_t dq_expected(fixture.state->get_nv());
                pinocchio::difference(*fixture.model,
                                      x0.head(fixture.state->get_nq()),
                                      x1.head(fixture.state->get_nq()),
                                      dq_expected);

                for (int i = 0; i < fixture.state->get_nv(); ++i)
                {
                    REQUIRE_THAT(dx(i), WithinAbs(dq_expected(i), tolerance<double>()));
                }

                // Verify velocity difference
                typename RS::VectorNv_t v0 = x0.tail(fixture.state->get_nv());
                typename RS::VectorNv_t v1 = x1.tail(fixture.state->get_nv());
                typename RS::VectorNv_t dv = dx.tail(fixture.state->get_nv());

                for (int i = 0; i < fixture.state->get_nv(); ++i)
                {
                    REQUIRE_THAT(dv(i), WithinAbs(v1(i) - v0(i), tolerance<double>()));
                }
            }

            SECTION("Integrate operation")
            {
                typename RS::VectorNdx_t dx = fixture.random_velocity();
                typename RS::VectorNx_t x_out(fixture.state->get_nx());
                fixture.state->integrate(x0, dx, x_out);

                // Verify configuration integration
                typename RS::VectorNq_t q_expected(fixture.state->get_nq());
                pinocchio::integrate(*fixture.model,
                                     x0.head(fixture.state->get_nq()),
                                     dx.head(fixture.state->get_nv()),
                                     q_expected);

                for (int i = 0; i < fixture.state->get_nq(); ++i)
                {
                    REQUIRE_THAT(x_out(i), WithinAbs(q_expected(i), tolerance<double>()));
                }

                // Verify velocity integration
                typename RS::VectorNv_t v0 = x0.tail(fixture.state->get_nv());
                typename RS::VectorNv_t dv = dx.tail(fixture.state->get_nv());
                typename RS::VectorNv_t v_out = x_out.tail(fixture.state->get_nv());

                for (int i = 0; i < fixture.state->get_nv(); ++i)
                {
                    REQUIRE_THAT(v_out(i), WithinAbs(v0(i) + dv(i), tolerance<double>()));
                }
            }

            SECTION("Round-trip consistency")
            {
                // Test: integrate(x0, diff(x0, x1)) ≈ x1
                typename RS::VectorNdx_t dx = fixture.state->diff_dx(x0, x1);
                typename RS::VectorNx_t x1_reconstructed = fixture.state->integrate_x(x0, dx);

                // Need special handling for quaternion comparison
                typename RS::VectorNq_t q1 = x1.head(fixture.state->get_nq());
                typename RS::VectorNq_t q1_recon = x1_reconstructed.head(fixture.state->get_nq());

                REQUIRE(pinocchio::isSameConfiguration(*fixture.model, q1, q1_recon, tolerance<double>() * 100));

                // Velocity comparison
                for (int i = fixture.state->get_nq(); i < fixture.state->get_nx(); ++i)
                {
                    REQUIRE_THAT(x1_reconstructed(i), WithinAbs(x1(i), tolerance<double>() * 10));
                }
            }
        }
    }
}

// Test Jacobian computations
TEST_CASE("StateMultibodyTpl - Jacobian Operations", "[multibody][jacobians]")
{
    using BS = BasicSpecTpl<double, double, 0>;
    using RS = RobotSpecTpl<BS, 7, 12, 6, 12, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

    MultibodyStateTestFixture<RS> fixture;
    fixture.setup_robot({"huron", 12, 12, true});

    SECTION("Jdiff correctness")
    {
        typename RS::VectorNx_t x0 = fixture.random_state();
        typename RS::VectorNx_t x1 = fixture.random_state();

        typename RS::MatrixNdx_t Jfirst(fixture.state->get_ndx(), fixture.state->get_ndx());
        typename RS::MatrixNdx_t Jsecond(fixture.state->get_ndx(), fixture.state->get_ndx());

        // Test both mode
        fixture.state->Jdiff(x0, x1, Jfirst, Jsecond, both);

        // Test individual modes
        SECTION("First component only")
        {
            typename RS::MatrixNdx_t J1(fixture.state->get_ndx(), fixture.state->get_ndx());
            typename RS::MatrixNdx_t J2(fixture.state->get_ndx(), fixture.state->get_ndx());
            J2.setZero();

            fixture.state->Jdiff(x0, x1, J1, J2, first);

            // J1 should match Jfirst, J2 should be zero
            for (int i = 0; i < J1.rows(); ++i)
            {
                for (int j = 0; j < J1.cols(); ++j)
                {
                    REQUIRE_THAT(J1(i, j), WithinAbs(Jfirst(i, j), tolerance<double>()));
                    REQUIRE_THAT(J2(i, j), WithinAbs(0.0, tolerance<double>()));
                }
            }
        }
    }

    SECTION("Jintegrate correctness")
    {
        typename RS::VectorNx_t x = fixture.random_state();
        typename RS::VectorNdx_t dx = fixture.random_velocity();

        typename RS::MatrixNdx_t Jfirst(fixture.state->get_ndx(), fixture.state->get_ndx());
        typename RS::MatrixNdx_t Jsecond(fixture.state->get_ndx(), fixture.state->get_ndx());

        // Test setto mode
        fixture.state->Jintegrate(x, dx, Jfirst, Jsecond, both, setto);

        // Test addto mode
        SECTION("Addto operation")
        {
            typename RS::MatrixNdx_t J1_add = RS::MatrixNdx_t::Ones(fixture.state->get_ndx(), fixture.state->get_ndx());
            typename RS::MatrixNdx_t J2_add = RS::MatrixNdx_t::Ones(fixture.state->get_ndx(), fixture.state->get_ndx());

            fixture.state->Jintegrate(x, dx, J1_add, J2_add, both, addto);

            // Should be original + Jacobian
            for (int i = 0; i < J1_add.rows(); ++i)
            {
                for (int j = 0; j < J1_add.cols(); ++j)
                {
                    REQUIRE_THAT(J1_add(i, j), WithinAbs(1.0 + Jfirst(i, j), tolerance<double>()));
                    REQUIRE_THAT(J2_add(i, j), WithinAbs(1.0 + Jsecond(i, j), tolerance<double>()));
                }
            }
        }

        // Test rmfrom mode
        SECTION("Rmfrom operation")
        {
            typename RS::MatrixNdx_t J1_rm = RS::MatrixNdx_t::Ones(fixture.state->get_ndx(), fixture.state->get_ndx());
            typename RS::MatrixNdx_t J2_rm = RS::MatrixNdx_t::Ones(fixture.state->get_ndx(), fixture.state->get_ndx());

            fixture.state->Jintegrate(x, dx, J1_rm, J2_rm, both, rmfrom);

            // Should be original - Jacobian
            for (int i = 0; i < J1_rm.rows(); ++i)
            {
                for (int j = 0; j < J1_rm.cols(); ++j)
                {
                    REQUIRE_THAT(J1_rm(i, j), WithinAbs(1.0 - Jfirst(i, j), tolerance<double>() * 10));
                    REQUIRE_THAT(J2_rm(i, j), WithinAbs(1.0 - Jsecond(i, j), tolerance<double>() * 10));
                }
            }
        }
    }

    SECTION("JintegrateTransport")
    {
        typename RS::VectorNx_t x = fixture.random_state();
        typename RS::VectorNdx_t dx = fixture.random_velocity() * 0.1; // Small displacement

        // Create a random matrix to transport
        typename RS::MatrixNv6_t Jin = RS::MatrixNv6_t::Random(fixture.state->get_nv(), 6);
        typename RS::MatrixNv6_t Jin_copy = Jin;

        // Transport for first argument
        fixture.state->JintegrateTransport(x, dx, Jin, first);

        // Verify transport operation via Pinocchio
        typename RS::MatrixNv6_t Jin_expected = Jin_copy;
        pinocchio::dIntegrateTransport(*fixture.model,
                                       x.head(fixture.state->get_nq()),
                                       dx.head(fixture.state->get_nv()),
                                       Jin_expected,
                                       pinocchio::ARG0);

        REQUIRE_THAT((Jin - Jin_expected).norm(), WithinAbs(0.0, tolerance<double>() * 100));
    }
}

// Test bounds management
TEST_CASE("StateMultibodyTpl - Bounds Management", "[multibody][bounds]")
{
    using BS = BasicSpecTpl<double, double, 0>;
    using RS = RobotSpecTpl<BS, 7, 30, 6, 30, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

    MultibodyStateTestFixture<RS> fixture;
    fixture.setup_robot({"atlas", 30, 30, true});

    SECTION("Initial bounds from model")
    {
        typename RS::VectorNx_t lb = fixture.state->get_lb();
        typename RS::VectorNx_t ub = fixture.state->get_ub();

        REQUIRE(lb.size() == fixture.state->get_nx());
        REQUIRE(ub.size() == fixture.state->get_nx());

        // Floating base bounds should be max, NOT infinity
        for (int i = 0; i < fixture.state->get_nqb(); ++i)
        {
            REQUIRE(lb(i) == -std::numeric_limits<double>::max());
            REQUIRE(ub(i) == std::numeric_limits<double>::max());
        }

        // Joint bounds should match model
        typename RS::VectorNq_t model_lb = fixture.model->lowerPositionLimit;
        typename RS::VectorNq_t model_ub = fixture.model->upperPositionLimit;

        for (int i = 0; i < fixture.state->get_nqj(); ++i)
        {
            int model_idx = fixture.state->get_nqb() + i;
            REQUIRE_THAT(lb(fixture.state->get_nqb() + i), WithinAbs(model_lb(model_idx), tolerance<double>()));
            REQUIRE_THAT(ub(fixture.state->get_nqb() + i), WithinAbs(model_ub(model_idx), tolerance<double>()));
        }

        // Velocity bounds
        typename RS::VectorNv_t vel_lb = lb.tail(fixture.state->get_nv());
        typename RS::VectorNv_t vel_ub = ub.tail(fixture.state->get_nv());

        for (int i = 0; i < fixture.state->get_nv(); ++i)
        {
            REQUIRE_THAT(vel_lb(i), WithinAbs(-fixture.model->velocityLimit(i), tolerance<double>()));
            REQUIRE_THAT(vel_ub(i), WithinAbs(fixture.model->velocityLimit(i), tolerance<double>()));
        }
    }

    SECTION("Bounds modification")
    {
        typename RS::VectorNx_t new_lb = RS::VectorNx_t::Constant(fixture.state->get_nx(), -1.0);
        typename RS::VectorNx_t new_ub = RS::VectorNx_t::Constant(fixture.state->get_nx(), 1.0);

        fixture.state->set_lb(new_lb);
        fixture.state->set_ub(new_ub);

        typename RS::VectorNx_t lb = fixture.state->get_lb();
        typename RS::VectorNx_t ub = fixture.state->get_ub();

        for (int i = 0; i < fixture.state->get_nx(); ++i)
        {
            REQUIRE_THAT(lb(i), WithinAbs(-1.0, tolerance<double>()));
            REQUIRE_THAT(ub(i), WithinAbs(1.0, tolerance<double>()));
        }
    }
}

// Test without floating base
TEST_CASE("StateMultibodyTpl - Fixed Base Robots", "[multibody][fixed_base]")
{
    using BS = BasicSpecTpl<double, double, 0>;
    using RS = RobotSpecTpl<BS, 0, 12, 0, 12, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

    MultibodyStateTestFixture<RS> fixture;

    // Load go1 without floating base
    fixture.model = std::make_unique<typename RS::RobotModel_t>();
    pinocchio::urdf::buildModel(robot_specs[2].urdf_path(), *fixture.model);
    fixture.data = std::make_unique<typename RS::RobotData_t>(*fixture.model);
    fixture.state = std::make_unique<typename RS::State_t>(*fixture.model);

    REQUIRE(fixture.state->get_nqb() == 0);
    REQUIRE(fixture.state->get_nvb() == 0);
    REQUIRE(fixture.state->get_nq() == fixture.state->get_nqj());
    REQUIRE(fixture.state->get_nv() == fixture.state->get_nvj());

    // Test basic operations still work
    typename RS::VectorNx_t x0 = fixture.state->zero();
    typename RS::VectorNx_t x1 = fixture.state->rand();
    typename RS::VectorNdx_t dx = fixture.state->diff_dx(x0, x1);
    typename RS::VectorNx_t x1_recon = fixture.state->integrate_x(x0, dx);

    // Should reconstruct x1
    for (int i = 0; i < fixture.state->get_nx(); ++i)
    {
        REQUIRE_THAT(x1_recon(i), WithinAbs(x1(i), tolerance<double>() * 100));
    }
}

// Template specialization tests - only test with double for now due to Pinocchio URDF parser limitations
TEST_CASE("StateMultibodyTpl - Double Scalar Type", "[multibody][scalar_types]")
{
    using Scalar = double;
    using BS = BasicSpecTpl<Scalar, Scalar, 0>;
    using RS = RobotSpecTpl<BS, 7, 12, 6, 12, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

    MultibodyStateTestFixture<RS> fixture;

    // Build model with appropriate scalar type
    fixture.model = std::make_unique<typename RS::RobotModel_t>();
    pinocchio::urdf::buildModel(robot_specs[1].urdf_path(), pinocchio::JointModelFreeFlyer(), *fixture.model);
    fixture.data = std::make_unique<typename RS::RobotData_t>(*fixture.model);
    fixture.state = std::make_unique<typename RS::State_t>(*fixture.model);

    SECTION("Basic operations with double")
    {
        typename RS::VectorNx_t x0 = fixture.state->zero();
        typename RS::VectorNx_t x1 = fixture.state->rand();

        // Test diff/integrate cycle
        typename RS::VectorNdx_t dx(fixture.state->get_ndx());
        fixture.state->diff(x0, x1, dx);

        typename RS::VectorNx_t x_out(fixture.state->get_nx());
        fixture.state->integrate(x0, dx, x_out);

        // Check reconstruction with appropriate tolerance
        typename RS::VectorNq_t q1 = x1.head(fixture.state->get_nq());
        typename RS::VectorNq_t q_out = x_out.head(fixture.state->get_nq());

        REQUIRE(pinocchio::isSameConfiguration(*fixture.model, q1, q_out, tolerance<Scalar>() * 1000));
    }
}

// Edge cases and error handling
TEST_CASE("StateMultibodyTpl - Edge Cases", "[multibody][edge_cases]")
{
    using BS = BasicSpecTpl<double, double, 0>;
    using RS = RobotSpecTpl<BS, 7, 12, 6, 12, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

    MultibodyStateTestFixture<RS> fixture;
    fixture.setup_robot({"go1", 12, 12, true});

    SECTION("Zero velocity integration")
    {
        typename RS::VectorNx_t x = fixture.random_state();
        typename RS::VectorNdx_t dx_zero = RS::VectorNdx_t::Zero(fixture.state->get_ndx());

        typename RS::VectorNx_t x_out = fixture.state->integrate_x(x, dx_zero);

        // Should return same state
        for (int i = 0; i < fixture.state->get_nx(); ++i)
        {
            REQUIRE_THAT(x_out(i), WithinAbs(x(i), tolerance<double>()));
        }
    }

    SECTION("Identity diff")
    {
        typename RS::VectorNx_t x = fixture.random_state();
        typename RS::VectorNdx_t dx = fixture.state->diff_dx(x, x);

        // Should be zero
        for (int i = 0; i < fixture.state->get_ndx(); ++i)
        {
            REQUIRE_THAT(dx(i), WithinAbs(0.0, tolerance<double>() * 10));
        }
    }

    SECTION("Large displacements")
    {
        typename RS::VectorNx_t x0 = fixture.state->zero();
        typename RS::VectorNdx_t dx_large = fixture.random_velocity() * 10.0;

        // Should still produce valid configuration
        typename RS::VectorNx_t x1 = fixture.state->integrate_x(x0, dx_large);
        typename RS::VectorNq_t q1 = x1.head(fixture.state->get_nq());

        REQUIRE(pinocchio::isNormalized(*fixture.model, q1, tolerance<double>() * 100));
    }
}

// Convenience method tests
TEST_CASE("StateMultibodyTpl - Convenience Methods", "[multibody][convenience]")
{
    using BS = BasicSpecTpl<double, double, 0>;
    using RS = RobotSpecTpl<BS, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

    MultibodyStateTestFixture<RS> fixture;
    fixture.setup_robot({"go1", 12, 12, true});

    SECTION("Jdiff_Js helper")
    {
        typename RS::VectorNx_t x0 = fixture.random_state();
        typename RS::VectorNx_t x1 = fixture.random_state();

        // Test both mode
        std::vector<typename RS::MatrixNdx_t> jacs_both = fixture.state->Jdiff_Js(x0, x1, both);
        REQUIRE(jacs_both.size() == 2);

        // Test first mode
        std::vector<typename RS::MatrixNdx_t> jacs_first = fixture.state->Jdiff_Js(x0, x1, first);
        REQUIRE(jacs_first.size() == 1);

        // Test second mode
        std::vector<typename RS::MatrixNdx_t> jacs_second = fixture.state->Jdiff_Js(x0, x1, second);
        REQUIRE(jacs_second.size() == 1);

        // Verify consistency
        typename RS::MatrixNdx_t J1(fixture.state->get_ndx(), fixture.state->get_ndx());
        typename RS::MatrixNdx_t J2(fixture.state->get_ndx(), fixture.state->get_ndx());
        fixture.state->Jdiff(x0, x1, J1, J2, both);

        for (int i = 0; i < J1.rows(); ++i)
        {
            for (int j = 0; j < J1.cols(); ++j)
            {
                REQUIRE_THAT(jacs_both[0](i, j), WithinAbs(J1(i, j), tolerance<double>()));
                REQUIRE_THAT(jacs_both[1](i, j), WithinAbs(J2(i, j), tolerance<double>()));
            }
        }
    }

    SECTION("Jintegrate_Js helper")
    {
        typename RS::VectorNx_t x = fixture.random_state();
        typename RS::VectorNdx_t dx = fixture.random_velocity();

        std::vector<typename RS::MatrixNdx_t> jacs = fixture.state->Jintegrate_Js(x, dx, both);
        REQUIRE(jacs.size() == 2);

        // Verify against direct call
        typename RS::MatrixNdx_t J1(fixture.state->get_ndx(), fixture.state->get_ndx());
        typename RS::MatrixNdx_t J2(fixture.state->get_ndx(), fixture.state->get_ndx());
        fixture.state->Jintegrate(x, dx, J1, J2, both, setto);

        for (int i = 0; i < J1.rows(); ++i)
        {
            for (int j = 0; j < J1.cols(); ++j)
            {
                REQUIRE_THAT(jacs[0](i, j), WithinAbs(J1(i, j), tolerance<double>()));
                REQUIRE_THAT(jacs[1](i, j), WithinAbs(J2(i, j), tolerance<double>()));
            }
        }
    }
}

// Performance characteristics test (compile-time vs runtime)
TEST_CASE("StateMultibodyTpl - Performance Characteristics", "[multibody][performance]")
{
    SECTION("Compile-time optimization verification")
    {
        using BS_fixed = BasicSpecTpl<double, double, 0>;
        using RS_fixed = RobotSpecTpl<BS_fixed, 7, 12, 6, 12, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

        using BS_dynamic = BasicSpecTpl<double, double, 0>;
        using RS_dynamic = RobotSpecTpl<BS_dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

        // Verify compile-time sizes
        static_assert(RS_fixed::NQ == 19);
        static_assert(RS_fixed::NV == 18);
        static_assert(RS_fixed::NX == 37);
        static_assert(RS_fixed::NDX == 36);

        // Verify dynamic sizes must be set at runtime
        static_assert(RS_dynamic::NQ == detail::Dynamic);
        static_assert(RS_dynamic::NV == detail::Dynamic);
    }
}
