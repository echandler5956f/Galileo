#ifndef __galileo_tests_helpers_catch_eigen_matchers_hpp__
#define __galileo_tests_helpers_catch_eigen_matchers_hpp__

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <Eigen/Dense>
#include <string>
#include <sstream>

namespace GalileoMatchers {

template <typename EigenType>
class EigenApproxMatcher : public Catch::Matchers::MatcherBase<EigenType> {
    EigenType m_expected;
    double m_epsilon;

public:
    EigenApproxMatcher(const EigenType& expected, double epsilon)
        : m_expected(expected), m_epsilon(epsilon) {}

    bool match(const EigenType& actual) const override {
        if (m_expected.rows() != actual.rows() || m_expected.cols() != actual.cols()) {
            return false;
        }
        return actual.isApprox(m_expected, m_epsilon);
    }

    std::string describe() const override {
        std::stringstream ss;
        ss << "is approximately equal to:\n" << m_expected << "\nwith tolerance " << m_epsilon;
        return ss.str();
    }
};

template <typename EigenType>
inline EigenApproxMatcher<EigenType> Approx(const EigenType& expected, double epsilon = 1e-8) {
    return EigenApproxMatcher<EigenType>(expected, epsilon);
}

} // namespace GalileoMatchers

#endif // __galileo_tests_helpers_catch_eigen_matchers_hpp__ 