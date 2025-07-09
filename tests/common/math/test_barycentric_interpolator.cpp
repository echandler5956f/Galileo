#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/common/math/barycentric-interpolator.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <sstream>
#include <chrono>

using namespace galileo;
using namespace Catch::Matchers;

template<typename NumScalar>
constexpr NumScalar TOLERANCE = std::numeric_limits<NumScalar>::epsilon() * 100;
