#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "galileo/common/math/radau-IIA.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace galileo::math;
using namespace Catch::Matchers;

TEST_CASE("RadauIIATpl", "[polynomial]")
{
    RadauIIATpl<double, 6, Eigen::ColMajor> poly;
    poly.compute_terms();

    Eigen::IOFormat CleanFmt(16, 0, ", ", "\n", "[", "]");
    std::cout << "A: \n" << poly.get_coefficients().format(CleanFmt) << std::endl;
    std::cout << "b: \n" << poly.get_weights().format(CleanFmt) << std::endl;
    std::cout << "c: \n" << poly.get_nodes().format(CleanFmt) << std::endl;
    std::cout << "W: \n" << poly.get_W().format(CleanFmt) << std::endl;
    std::cout << "X: \n" << poly.get_X().format(CleanFmt) << std::endl;
}
