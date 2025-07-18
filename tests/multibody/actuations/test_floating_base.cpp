#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/multibody/actuations/implementations/actuation-floating-base.hpp"
#include "galileo/core/states/implementations/state-euclidean.hpp"
#include "galileo/core/basic-spec.hpp"
#include "galileo/multibody/robot-spec.hpp"
#include "galileo/common/meta/dimension.hpp"
#include "galileo/common/meta/eigen.hpp"

#include <Eigen/Dense>
#include <memory>
#include <limits>
#include <random>
#include <iostream>

using namespace galileo;
using namespace Catch::Matchers;

template <typename Scalar>
constexpr Scalar TOLERANCE = std::numeric_limits<Scalar>::epsilon() * 1000;

// Test fixture for common test setup
template <typename Scalar>
class FloatingBaseTestFixture
{
public:
    using BS = BasicSpecTpl<Scalar, Scalar, 0>;

    // Static configuration for testing
    static constexpr int NQb = 6;
    static constexpr int NQj = 12;
    static constexpr int NVb = 6;
    static constexpr int NVj = 12;
    static constexpr int NRotors = 0;

    using RS = RobotSpecTpl<BS, NQb, NQj, NVb, NVj, NRotors, StateEuclideanTpl, ActuationFloatingBaseTpl>;
    using ActuationModel_t = ActuationModelFloatingBaseTpl<RS>;
    using ActuationData_t = ActuationDataTpl<RS>;

    FloatingBaseTestFixture()
    {
        // Set up robot spec with bounds
        typename RS::VectorNx_t lb = RS::VectorNx_t::Constant(RS::NX, Scalar(-10.0));
        typename RS::VectorNx_t ub = RS::VectorNx_t::Constant(RS::NX, Scalar(10.0));

        // Create shared state
        state_ = std::make_shared<StateEuclideanTpl<RS>>(rs_, lb, ub);

        // Create actuation model
        actuation_model_ = std::make_unique<ActuationModel_t>(state_);

        // Create test vectors
        x_ = RS::VectorNx_t::Random(RS::NX);
        u_ = RS::VectorNua_t::Random(RS::NUa);
        tau_ = RS::VectorNv_t::Random(RS::NV);
    }

    RS rs_;
    std::shared_ptr<StateEuclideanTpl<RS>> state_;
    std::unique_ptr<ActuationModel_t> actuation_model_;

    typename RS::VectorNx_t x_;
    typename RS::VectorNua_t u_;
    typename RS::VectorNv_t tau_;
};

// Helper function for testing matrix properties
template <typename MatrixType>
void test_matrix_properties(const MatrixType& matrix, int expected_rows, int expected_cols, bool should_be_zero = false)
{
    REQUIRE(matrix.rows() == expected_rows);
    REQUIRE(matrix.cols() == expected_cols);

    if (should_be_zero) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                REQUIRE_THAT(static_cast<double>(matrix(i, j)), WithinAbs(0.0, TOLERANCE<double>));
            }
        }
    }
}

struct TestDataCollector
{
};

// Test basic construction and properties
TEST_CASE("ActuationModelFloatingBaseTpl - Basic Construction", "[floating_base][construction]")
{
    SECTION("Fixed dimensions construction")
    {
        FloatingBaseTestFixture<double> fixture;

        // Test that the model was constructed successfully
        REQUIRE(fixture.actuation_model_ != nullptr);

        // Test dimension accessors
        REQUIRE(fixture.actuation_model_->get_state()->get_nua() == fixture.rs_.NUa);
        REQUIRE(fixture.actuation_model_->get_state()->get_nua() == FloatingBaseTestFixture<double>::NVj + FloatingBaseTestFixture<double>::NRotors);
        REQUIRE(fixture.actuation_model_->get_state()->get_nua() == 12); // 12 + 0

        // Test state access
        REQUIRE(fixture.actuation_model_->get_state() == fixture.state_);
    }

    SECTION("Different floating base configurations")
    {
        using BS = BasicSpecTpl<double, double, 0>;

        // Test with no rotors
        {
            using RS_NoRotors = RobotSpecTpl<BS, 6, 12, 6, 12, 0, StateEuclideanTpl, ActuationFloatingBaseTpl>;
            RS_NoRotors rs;
            typename RS_NoRotors::VectorNx_t lb = RS_NoRotors::VectorNx_t::Constant(RS_NoRotors::NX, -5.0);
            typename RS_NoRotors::VectorNx_t ub = RS_NoRotors::VectorNx_t::Constant(RS_NoRotors::NX, 5.0);

            auto state = std::make_shared<StateEuclideanTpl<RS_NoRotors>>(rs, lb, ub);
            ActuationModelFloatingBaseTpl<RS_NoRotors> model(state);

            REQUIRE(model.get_state()->get_nua() == 12); // Only joint actuations
        }

        // Test with many rotors
        {
            using RS_ManyRotors = RobotSpecTpl<BS, 6, 12, 6, 12, 4, StateEuclideanTpl, ActuationFloatingBaseTpl>;
            RS_ManyRotors rs;
            typename RS_ManyRotors::VectorNx_t lb = RS_ManyRotors::VectorNx_t::Constant(RS_ManyRotors::NX, -5.0);
            typename RS_ManyRotors::VectorNx_t ub = RS_ManyRotors::VectorNx_t::Constant(RS_ManyRotors::NX, 5.0);

            auto state = std::make_shared<StateEuclideanTpl<RS_ManyRotors>>(rs, lb, ub);
            ActuationModelFloatingBaseTpl<RS_ManyRotors> model(state);

            REQUIRE(model.get_state()->get_nua() == 16); // 12 joint + 4 rotor actuations
        }
    }
}

// Test createData method
TEST_CASE("ActuationModelFloatingBaseTpl - CreateData", "[floating_base][create_data]")
{
    SECTION("Fixed dimensions data creation")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;

        // Create data
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Test data dimensions
        REQUIRE(data.tau.size() == fixture.rs_.NV);
        REQUIRE(data.u.size() == fixture.rs_.NUa);
        REQUIRE(data.dtau_dx.rows() == fixture.rs_.NV);
        REQUIRE(data.dtau_dx.cols() == fixture.rs_.NDX);
        REQUIRE(data.dtau_du.rows() == fixture.rs_.NV);
        REQUIRE(data.dtau_du.cols() == fixture.rs_.NUa);
        REQUIRE(data.Mtau.rows() == fixture.rs_.NUa);
        REQUIRE(data.Mtau.cols() == fixture.rs_.NV);
        REQUIRE(data.tau_set.size() == fixture.rs_.NV);

        // Test initial values
        for (int i = 0; i < data.tau.size(); ++i) {
            REQUIRE_THAT(data.tau(i), WithinAbs(0.0, TOLERANCE<double>));
        }

        for (int i = 0; i < data.u.size(); ++i) {
            REQUIRE_THAT(data.u(i), WithinAbs(0.0, TOLERANCE<double>));
        }

        // Test that dtau_dx is zero
        test_matrix_properties(data.dtau_dx, fixture.rs_.NV, fixture.rs_.NDX, true);

        // Test dtau_du structure - should have identity in the actuated part
        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            for (int j = 0; j < fixture.rs_.NUa; ++j) {
                if (i == j) {
                    REQUIRE_THAT(data.dtau_du(fixture.rs_.NVb + i, j), WithinAbs(1.0, TOLERANCE<double>));
                } else {
                    REQUIRE_THAT(data.dtau_du(fixture.rs_.NVb + i, j), WithinAbs(0.0, TOLERANCE<double>));
                }
            }
        }

        // Test Mtau structure - should have identity in the actuated part
        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            for (int j = 0; j < fixture.rs_.NV; ++j) {
                if (j == fixture.rs_.NVb + i) {
                    REQUIRE_THAT(data.Mtau(i, j), WithinAbs(1.0, TOLERANCE<double>));
                } else {
                    REQUIRE_THAT(data.Mtau(i, j), WithinAbs(0.0, TOLERANCE<double>));
                }
            }
        }

        // Test tau_set - floating base should be false, actuated should be true
        for (int i = 0; i < fixture.rs_.NVb; ++i) {
            REQUIRE_FALSE(data.tau_set(i));
        }
        for (int i = fixture.rs_.NVb; i < fixture.rs_.NV; ++i) {
            REQUIRE(data.tau_set(i));
        }
    }
}

// Test calc method
TEST_CASE("ActuationModelFloatingBaseTpl - Calc Method", "[floating_base][calc]")
{
    SECTION("Basic calc functionality")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Call calc
        fixture.actuation_model_->calc(data, fixture.x_, fixture.u_);

        // Test that tau is set correctly - tail should equal u
        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            REQUIRE_THAT(data.tau(fixture.rs_.NVb + i), WithinAbs(fixture.u_(i), TOLERANCE<double>));
        }

        // Test that floating base part remains zero
        for (int i = 0; i < fixture.rs_.NVb; ++i) {
            REQUIRE_THAT(data.tau(i), WithinAbs(0.0, TOLERANCE<double>));
        }
    }

    SECTION("Calc with different vector types")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Test with dynamic vectors
        Eigen::VectorXd x_dynamic = fixture.x_;
        Eigen::VectorXd u_dynamic = fixture.u_;

        fixture.actuation_model_->calc(data, x_dynamic, u_dynamic);

        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            REQUIRE_THAT(data.tau(fixture.rs_.NVb + i), WithinAbs(fixture.u_(i), TOLERANCE<double>));
        }
    }

    SECTION("Calc with zero input")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        typename FloatingBaseTestFixture<double>::RS::VectorNua_t u_zero =
            FloatingBaseTestFixture<double>::RS::VectorNua_t::Zero(fixture.rs_.NUa);

        fixture.actuation_model_->calc(data, fixture.x_, u_zero);

        // All tau should be zero
        for (int i = 0; i < data.tau.size(); ++i) {
            REQUIRE_THAT(data.tau(i), WithinAbs(0.0, TOLERANCE<double>));
        }
    }
}

// Test calcDiff method
TEST_CASE("ActuationModelFloatingBaseTpl - CalcDiff Method", "[floating_base][calc_diff]")
{
    SECTION("CalcDiff does not modify jacobians")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Store original jacobian values
        auto dtau_dx_orig = data.dtau_dx;
        auto dtau_du_orig = data.dtau_du;
        auto Mtau_orig = data.Mtau;

        fixture.actuation_model_->calcDiff(data, fixture.x_, fixture.u_);

        // Jacobians should remain unchanged (floating base has constant jacobians)
        for (int i = 0; i < data.dtau_dx.rows(); ++i) {
            for (int j = 0; j < data.dtau_dx.cols(); ++j) {
                REQUIRE_THAT(data.dtau_dx(i, j), WithinAbs(dtau_dx_orig(i, j), TOLERANCE<double>));
            }
        }

        for (int i = 0; i < data.dtau_du.rows(); ++i) {
            for (int j = 0; j < data.dtau_du.cols(); ++j) {
                REQUIRE_THAT(data.dtau_du(i, j), WithinAbs(dtau_du_orig(i, j), TOLERANCE<double>));
            }
        }

        for (int i = 0; i < data.Mtau.rows(); ++i) {
            for (int j = 0; j < data.Mtau.cols(); ++j) {
                REQUIRE_THAT(data.Mtau(i, j), WithinAbs(Mtau_orig(i, j), TOLERANCE<double>));
            }
        }
    }
}

// Test commands method
TEST_CASE("ActuationModelFloatingBaseTpl - Commands Method", "[floating_base][commands]")
{
    SECTION("Commands extraction from tau")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Set tau values
        for (int i = 0; i < fixture.rs_.NVb; ++i) {
            fixture.tau_(i) = 100.0 + i; // Should be ignored
        }
        for (int i = fixture.rs_.NVb; i < fixture.rs_.NV; ++i) {
            fixture.tau_(i) = 10.0 + i; // Should be extracted
        }

        fixture.actuation_model_->commands(data, fixture.x_, fixture.tau_);

        // Check that u contains the tail of tau
        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            REQUIRE_THAT(data.u(i), WithinAbs(fixture.tau_(fixture.rs_.NVb + i), TOLERANCE<double>));
        }
    }

    SECTION("Commands with different tau configurations")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Test with all zeros
        typename FloatingBaseTestFixture<double>::RS::VectorNv_t tau_zero =
            FloatingBaseTestFixture<double>::RS::VectorNv_t::Zero(fixture.rs_.NV);

        fixture.actuation_model_->commands(data, fixture.x_, tau_zero);

        for (int i = 0; i < data.u.size(); ++i) {
            REQUIRE_THAT(data.u(i), WithinAbs(0.0, TOLERANCE<double>));
        }

        // Test with specific pattern
        typename FloatingBaseTestFixture<double>::RS::VectorNv_t tau_pattern =
            FloatingBaseTestFixture<double>::RS::VectorNv_t::Zero(fixture.rs_.NV);

        for (int i = fixture.rs_.NVb; i < fixture.rs_.NV; ++i) {
            tau_pattern(i) = i * 2.5;
        }

        fixture.actuation_model_->commands(data, fixture.x_, tau_pattern);

        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            REQUIRE_THAT(data.u(i), WithinAbs((fixture.rs_.NVb + i) * 2.5, TOLERANCE<double>));
        }
    }
}

// Test torqueTransform method
TEST_CASE("ActuationModelFloatingBaseTpl - TorqueTransform Method", "[floating_base][torque_transform]")
{
    SECTION("TorqueTransform preserves jacobians")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Store original values
        auto Mtau_orig = data.Mtau;

        fixture.actuation_model_->torqueTransform(data, fixture.x_, fixture.u_);

        // Mtau should remain unchanged (floating base has constant transform)
        for (int i = 0; i < data.Mtau.rows(); ++i) {
            for (int j = 0; j < data.Mtau.cols(); ++j) {
                REQUIRE_THAT(data.Mtau(i, j), WithinAbs(Mtau_orig(i, j), TOLERANCE<double>));
            }
        }
    }
}

// Test dynamic dimensions
TEST_CASE("ActuationModelFloatingBaseTpl - Dynamic Dimensions", "[floating_base][dynamic]")
{
    SECTION("Dynamic robot spec")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS_Dynamic = RobotSpecTpl<BS, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, detail::Dynamic, StateEuclideanTpl, ActuationFloatingBaseTpl>;

        RS_Dynamic rs;

        // Setup dynamic dimensions
        rs.nqb_dim_.set_value(6);
        rs.nqj_dim_.set_value(15);
        rs.nvb_dim_.set_value(6);
        rs.nvj_dim_.set_value(15);
        rs.nrotors_dim_.set_value(3);
        rs.nq_dim_.set_value(21);
        rs.nv_dim_.set_value(21);
        rs.nx_dim_.set_value(42);
        rs.ndx_dim_.set_value(42);
        rs.nua_dim_.set_value(18); // 15 + 3

        typename RS_Dynamic::VectorNx_t lb(42);
        typename RS_Dynamic::VectorNx_t ub(42);
        lb.setConstant(-5.0);
        ub.setConstant(5.0);

        auto state = std::make_shared<StateEuclideanTpl<RS_Dynamic>>(rs, lb, ub);
        ActuationModelFloatingBaseTpl<RS_Dynamic> model(state);

        // Test dimensions
        REQUIRE(model.get_state()->get_nua() == 18);

        // Test data creation
        TestDataCollector data_collector;
        auto data = model.createData(&data_collector);
        REQUIRE(data.tau.size() == 21);
        REQUIRE(data.u.size() == 18);
        REQUIRE(data.dtau_du.rows() == 21);
        REQUIRE(data.dtau_du.cols() == 18);

        // Test computation
        typename RS_Dynamic::VectorNx_t x = RS_Dynamic::VectorNx_t::Random(42);
        typename RS_Dynamic::VectorNua_t u = RS_Dynamic::VectorNua_t::Random(18);

        model.calc(data, x, u);

        // Check that actuated part is set correctly
        int start_idx = data.tau.size() - u.size();  // NV - NUa
        for (int i = 0; i < 18; ++i) {
            REQUIRE_THAT(data.tau(start_idx + i), WithinAbs(u(i), TOLERANCE<double>));
        }
    }
}

// Test with different scalar types
TEMPLATE_TEST_CASE("ActuationModelFloatingBaseTpl - Different Scalar Types", "[floating_base][scalar_types]", float, double, long double)
{
    using Scalar = TestType;
    constexpr Scalar tolerance = TOLERANCE<Scalar>;

    SECTION("Basic operations with different scalar types")
    {
        FloatingBaseTestFixture<Scalar> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Test calc
        fixture.actuation_model_->calc(data, fixture.x_, fixture.u_);

        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            REQUIRE_THAT(static_cast<double>(data.tau(fixture.rs_.NVb + i)),
                        WithinAbs(static_cast<double>(fixture.u_(i)), tolerance));
        }

        // Test commands
        fixture.actuation_model_->commands(data, fixture.x_, fixture.tau_);

        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            REQUIRE_THAT(static_cast<double>(data.u(i)),
                        WithinAbs(static_cast<double>(fixture.tau_(fixture.rs_.NVb + i)), tolerance));
        }
    }
}

// Test integration with full computation pipeline
TEST_CASE("ActuationModelFloatingBaseTpl - Integration Tests", "[floating_base][integration]")
{
    SECTION("Full computation pipeline")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data = fixture.actuation_model_->createData(&data_collector);

        // Full pipeline: calc -> commands -> calc again
        fixture.actuation_model_->calc(data, fixture.x_, fixture.u_);

        // Extract tau for commands
        typename FloatingBaseTestFixture<double>::RS::VectorNv_t tau_extracted = data.tau;

        // Clear u and recompute via commands
        data.u.setZero();
        fixture.actuation_model_->commands(data, fixture.x_, tau_extracted);

        // u should now match original
        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            REQUIRE_THAT(data.u(i), WithinAbs(fixture.u_(i), TOLERANCE<double>));
        }
    }

    SECTION("Consistency across multiple calls")
    {
        FloatingBaseTestFixture<double> fixture;
        TestDataCollector data_collector;
        auto data1 = fixture.actuation_model_->createData(&data_collector);
        auto data2 = fixture.actuation_model_->createData(&data_collector);

        // Same inputs should give same outputs
        fixture.actuation_model_->calc(data1, fixture.x_, fixture.u_);
        fixture.actuation_model_->calc(data2, fixture.x_, fixture.u_);

        for (int i = 0; i < data1.tau.size(); ++i) {
            REQUIRE_THAT(data1.tau(i), WithinAbs(data2.tau(i), TOLERANCE<double>));
        }

        // Different inputs should give different outputs
        typename FloatingBaseTestFixture<double>::RS::VectorNua_t u_different =
            fixture.u_ + FloatingBaseTestFixture<double>::RS::VectorNua_t::Ones(fixture.rs_.NUa);

        fixture.actuation_model_->calc(data2, fixture.x_, u_different);

        bool any_different = false;
        for (int i = fixture.rs_.NVb; i < fixture.rs_.NV; ++i) {
            if (std::abs(data1.tau(i) - data2.tau(i)) > TOLERANCE<double>) {
                any_different = true;
                break;
            }
        }
        REQUIRE(any_different);
    }
}

// Test dimension calculations with DimensionTpl
TEST_CASE("ActuationModelFloatingBaseTpl - DimensionTpl Integration", "[floating_base][dimensions]")
{
    SECTION("Compile-time dimension verification")
    {
        FloatingBaseTestFixture<double> fixture;

        // Test that compile-time dimensions match runtime
        static_assert(FloatingBaseTestFixture<double>::RS::NUa == FloatingBaseTestFixture<double>::NVj + FloatingBaseTestFixture<double>::NRotors);
        static_assert(FloatingBaseTestFixture<double>::RS::NV == FloatingBaseTestFixture<double>::NVb + FloatingBaseTestFixture<double>::NVj);
        static_assert(FloatingBaseTestFixture<double>::RS::NUa == 12);
        static_assert(FloatingBaseTestFixture<double>::RS::NV == 18);

        REQUIRE(fixture.actuation_model_->get_state()->get_nua() == FloatingBaseTestFixture<double>::RS::NUa);
        REQUIRE(fixture.state_->get_nv() == FloatingBaseTestFixture<double>::RS::NV);
        REQUIRE(fixture.state_->get_nua() == FloatingBaseTestFixture<double>::RS::NUa);
    }

    SECTION("Mixed compile-time and dynamic dimensions")
    {
        using BS = BasicSpecTpl<double, double, 0>;
        using RS_Mixed = RobotSpecTpl<BS, 6, detail::Dynamic, 6, detail::Dynamic, 0, StateEuclideanTpl, ActuationFloatingBaseTpl>;

        RS_Mixed rs;

        // Set dynamic dimensions
        rs.nqj_dim_.set_value(12);
        rs.nvj_dim_.set_value(12);
        rs.nq_dim_.set_value(18);
        rs.nv_dim_.set_value(18);
        rs.nx_dim_.set_value(36);
        rs.ndx_dim_.set_value(36);
        rs.nua_dim_.set_value(12); // 12 + 0

        typename RS_Mixed::VectorNx_t lb(36);
        typename RS_Mixed::VectorNx_t ub(36);
        lb.setConstant(-5.0);
        ub.setConstant(5.0);

        auto state = std::make_shared<StateEuclideanTpl<RS_Mixed>>(rs, lb, ub);
        ActuationModelFloatingBaseTpl<RS_Mixed> model(state);

        // Test mixed dimensions work correctly
        REQUIRE(model.get_state()->get_nua() == 12);

        TestDataCollector data_collector;
        auto data = model.createData(&data_collector);

        // Test matrix dimensions
        REQUIRE(data.dtau_du.rows() == 18);
        REQUIRE(data.dtau_du.cols() == 12);
        REQUIRE(data.Mtau.rows() == 12);
        REQUIRE(data.Mtau.cols() == 18);

        // Test computation
        typename RS_Mixed::VectorNx_t x = RS_Mixed::VectorNx_t::Random(36);
        typename RS_Mixed::VectorNua_t u = RS_Mixed::VectorNua_t::Random(12);

        model.calc(data, x, u);

        for (int i = 0; i < 12; ++i) {
            int start_idx = data.tau.size() - u.size();  // NV - NUa
            REQUIRE_THAT(data.tau(start_idx + i), WithinAbs(u(i), TOLERANCE<double>));
        }
    }
}

// Test thread safety and const-correctness
TEST_CASE("ActuationModelFloatingBaseTpl - Thread Safety", "[floating_base][thread_safety]")
{
    SECTION("Const-correctness")
    {
        FloatingBaseTestFixture<double> fixture;
        const auto& const_model = *fixture.actuation_model_;

        // Test const methods
        REQUIRE(const_model.get_state()->get_nua() == fixture.rs_.NUa);

        // Test that const methods don't modify state
        TestDataCollector data_collector;
        auto data = const_model.createData(&data_collector);

        // These should all be const operations
        const_model.calc(data, fixture.x_, fixture.u_);
        const_model.calcDiff(data, fixture.x_, fixture.u_);
        const_model.commands(data, fixture.x_, fixture.tau_);
        const_model.torqueTransform(data, fixture.x_, fixture.u_);
    }

    SECTION("Multiple data objects")
    {
        FloatingBaseTestFixture<double> fixture;

        // Create multiple data objects
        TestDataCollector data_collector;
        auto data1 = fixture.actuation_model_->createData(&data_collector);
        auto data2 = fixture.actuation_model_->createData(&data_collector);

        // Different operations on different data should not interfere
        typename FloatingBaseTestFixture<double>::RS::VectorNua_t u1 = fixture.u_;
        typename FloatingBaseTestFixture<double>::RS::VectorNua_t u2 = fixture.u_ * 2.0;

        fixture.actuation_model_->calc(data1, fixture.x_, u1);
        fixture.actuation_model_->calc(data2, fixture.x_, u2);

        // Results should be independent
        for (int i = 0; i < fixture.rs_.NUa; ++i) {
            REQUIRE_THAT(data1.tau(fixture.rs_.NVb + i), WithinAbs(u1(i), TOLERANCE<double>));
            REQUIRE_THAT(data2.tau(fixture.rs_.NVb + i), WithinAbs(u2(i), TOLERANCE<double>));
        }
    }
}
