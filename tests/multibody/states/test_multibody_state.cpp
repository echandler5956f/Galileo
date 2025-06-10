#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include "galileo/multibody/states/multibody.hpp"
#include "galileo/multibody/actuations/floating-base.hpp"
#include "galileo/multibody/robot-spec.hpp"
#include "../../models/pendulum.hpp"
#include "../../helpers/catch_eigen_matchers.hpp"
#include "../../helpers/finite_difference.hpp"

using namespace galileo;
using namespace GalileoMatchers;

// Define a simple RobotSpec for the pendulum
// NQb=0, NQj=1, NVb=0, NVj=1, NRotors=0
using PendulumBasicSpec = BasicSpecTpl<double, double, 0>;
template <template <typename> class ActuationTpl>
using PendulumRobotSpecTpl = RobotSpecTpl<PendulumBasicSpec, 0, 1, 0, 1, 0, StateMultibodyTpl, ActuationTpl>;
using PendulumRobotSpec = PendulumRobotSpecTpl<ActuationModelFloatingBaseTpl>;

TEMPLATE_TEST_CASE("Multibody State Operations", "[multibody][state]", PendulumRobotSpec)
{
    using RobotSpec = TestType;
    using State = typename RobotSpec::State_t;
    using RobotModel = typename RobotSpec::RobotModel_t;
    using VectorNx = typename RobotSpec::VectorNx_t;
    using VectorNdx = typename RobotSpec::VectorNdx_t;
    using MatrixNdx = typename RobotSpec::MatrixNdx_t;
    using Scalar = typename RobotSpec::NumScalar;

    RobotModel model;
    tests::build_pendulum(model);
    State state(&model);

    VectorNx x0 = state.rand();
    VectorNx x1 = state.rand();
    VectorNdx dx = VectorNdx::Random(state.get_ndx());

    WHEN("diff and integrate are called")
    {
        VectorNdx dx_test(state.get_ndx());
        state.diff(x0, x1, dx_test);

        VectorNx x1_integrated(state.get_nx());
        state.integrate(x0, dx_test, x1_integrated);

        THEN("integrate(x0, diff(x0, x1)) should recover x1")
        {
            REQUIRE_THAT(x1_integrated, Approx(x1, 1e-9));
        }
    }

    WHEN("Jdiff is called")
    {
        MatrixNdx Jfirst_analytical(state.get_ndx(), state.get_ndx());
        MatrixNdx Jsecond_analytical(state.get_ndx(), state.get_ndx());
        state.Jdiff(x0, x1, Jfirst_analytical, Jsecond_analytical);

        THEN("the analytical Jacobians match the finite difference approximation")
        {
            MatrixNdx Jfirst_numerical(state.get_ndx(), state.get_ndx());
            auto func_J0 = [&](const VectorNx &x, VectorNdx &y_out)
            { state.diff(x, x1, y_out); };
            tests::compute_finite_difference_jacobian_manifold<State, decltype(func_J0), VectorNx, VectorNdx, MatrixNdx>(
                state, func_J0, x0, 1e-7, Jfirst_numerical);
            REQUIRE_THAT(Jfirst_analytical, Approx(Jfirst_numerical, 1e-5));

            MatrixNdx Jsecond_numerical(state.get_ndx(), state.get_ndx());
            auto func_J1 = [&](const VectorNx &x, VectorNdx &y_out)
            { state.diff(x0, x, y_out); };
            tests::compute_finite_difference_jacobian_manifold<State, decltype(func_J1), VectorNx, VectorNdx, MatrixNdx>(
                state, func_J1, x1, 1e-7, Jsecond_numerical);
            REQUIRE_THAT(Jsecond_analytical, Approx(Jsecond_numerical, 1e-5));
        }
    }

    WHEN("Jintegrate is called")
    {
        MatrixNdx Jfirst_analytical(state.get_ndx(), state.get_ndx());
        MatrixNdx Jsecond_analytical(state.get_ndx(), state.get_ndx());
        state.Jintegrate(x0, dx, Jfirst_analytical, Jsecond_analytical);

        THEN("the analytical Jacobians match the finite difference approximation")
        {
            MatrixNdx Jfirst_numerical(state.get_ndx(), state.get_ndx());
            VectorNx x_out(state.get_nx());
            VectorNdx diff_out(state.get_ndx());
            auto func_J0 = [&](const VectorNx &x, VectorNdx &y_out)
            {
                state.integrate(x, dx, x_out);
                state.diff(x0, x_out, y_out);
            };
            tests::compute_finite_difference_jacobian_manifold<State, decltype(func_J0), VectorNx, VectorNdx, MatrixNdx>(
                state, func_J0, x0, 1e-7, Jfirst_numerical);
            REQUIRE_THAT(Jfirst_analytical, Approx(Jfirst_numerical, 1e-5));

            MatrixNdx Jsecond_numerical(state.get_ndx(), state.get_ndx());
            auto func_J1 = [&](const VectorNdx &d, VectorNdx &y_out)
            {
                state.integrate(x0, d, x_out);
                state.diff(x0, x_out, y_out);
            };
            tests::compute_finite_difference_jacobian<decltype(func_J1), VectorNdx, VectorNdx, VectorNdx>(
                func_J1, dx, 1e-7, Jsecond_numerical);
            REQUIRE_THAT(Jsecond_analytical, Approx(Jsecond_numerical, 1e-5));
        }
    }
}