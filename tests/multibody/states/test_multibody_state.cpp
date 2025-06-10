#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include "galileo/multibody/states/multibody.hpp"
#include "galileo/multibody/actuations/floating-base.hpp"
#include "galileo/multibody/robot-spec.hpp"
#include "../../models/pendulum.hpp"
#include "../../helpers/catch_eigen_matchers.hpp"

using namespace galileo;
using namespace GalileoMatchers;

// Define a simple RobotSpec for the pendulum
// NQb=0, NQj=1, NVb=0, NVj=1, NRotors=0
using PendulumBasicSpec = BasicSpecTpl<double, double, 0>;
template <template<typename> class ActuationTpl>
using PendulumRobotSpecTpl = RobotSpecTpl<PendulumBasicSpec, 0, 1, 0, 1, 0, StateMultibodyTpl, ActuationTpl>;
using PendulumRobotSpec = PendulumRobotSpecTpl<ActuationModelFloatingBaseTpl>;


TEMPLATE_TEST_CASE("Multibody State Operations", "[multibody][state]", PendulumRobotSpec) {
    using RobotSpec = TestType;
    using State = typename RobotSpec::State_t;
    using RobotModel = typename RobotSpec::RobotModel_t;
    using VectorNx = typename RobotSpec::VectorNx_t;
    using VectorNdx = typename RobotSpec::VectorNdx_t;

    RobotModel model;
    tests::build_pendulum(model);
    State state(&model);

    GIVEN("Two random states x0 and x1") {
        VectorNx x0 = state.zero();
        VectorNx x1 = state.rand();
        
        WHEN("the state difference is computed and then integrated") {
            VectorNdx dx(state.get_ndx());
            state.diff(x0, x1, dx);

            VectorNx x1_integrated(state.get_nx());
            state.integrate(x0, dx, x1_integrated);

            THEN("the original state x1 is recovered") {
                REQUIRE_THAT(x1_integrated, Approx(x1, 1e-9));
            }
        }
    }

    // A more advanced test would be to check the Jacobians with finite differences.
    // This serves as a good starting point.
} 