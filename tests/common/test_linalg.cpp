#include <catch2/catch_test_macros.hpp>
#include "galileo/common/linalg/helpers.hpp"
#include "../helpers/catch_eigen_matchers.hpp"

using namespace galileo;
using namespace GalileoMatchers;

SCENARIO("Matrix properties can be checked and enforced", "[common][linalg]") {
    GIVEN("A symmetric matrix") {
        Eigen::Matrix3d M;
        M << 1, 2, 3,
             2, 5, 6,
             3, 6, 9;

        WHEN("isSymmetric is called") {
            THEN("it should return true") {
                REQUIRE(math::isSymmetric(M) == true);
            }
        }
    }

    GIVEN("A non-symmetric matrix") {
        Eigen::Matrix3d M;
        M << 1, 2, 3,
             4, 5, 6,
             7, 8, 9;

        WHEN("isSymmetric is called") {
            THEN("it should return false") {
                REQUIRE(math::isSymmetric(M) == false);
            }
        }

        WHEN("enforceSymmetric is called") {
            math::enforceSymmetric(M);
            THEN("the matrix becomes symmetric") {
                REQUIRE(math::isSymmetric(M) == true);

                Eigen::Matrix3d expected;
                expected << 1.0, 3.0, 5.0,
                            3.0, 5.0, 7.0,
                            5.0, 7.0, 9.0;
                REQUIRE_THAT(M, Approx(expected));
            }
        }
    }

    GIVEN("A positive definite matrix") {
        Eigen::Matrix2d M;
        M << 2, 1,
             1, 2;

        WHEN("isPositiveDefinite is called") {
            THEN("it should return true") {
                REQUIRE(math::isPositiveDefinite(M) == true);
            }
        }
    }

    GIVEN("A non-positive definite matrix") {
        Eigen::Matrix2d M;
        M << 1, 2,
             2, 1;

        WHEN("isPositiveDefinite is called") {
            THEN("it should return false") {
                REQUIRE(math::isPositiveDefinite(M) == false);
            }
        }

        WHEN("enforcePositiveDefinite is called") {
            math::enforcePositiveDefinite(M);
            THEN("the matrix becomes positive definite") {
                REQUIRE(math::isPositiveDefinite(M) == true);
            }
        }
    }
}
