#ifndef GALILEO_TESTING_UTILS_MATCHERS_EIGEN_MATCHERS_HPP
#define GALILEO_TESTING_UTILS_MATCHERS_EIGEN_MATCHERS_HPP

#include <catch2/matchers/catch_matchers_templated.hpp>
#include <Eigen/Core>
#include <sstream>

namespace galileo
{
    namespace testing
    {
        template <typename EigenType>
        struct EigenApproxMatcher : Catch::Matchers::MatcherBase<EigenType>
        {
            EigenApproxMatcher(const EigenType &expected, double tolerance)
                : m_expected(expected), m_tolerance(tolerance) {}

            bool match(const EigenType &actual) const override
            {
                if (m_expected.rows() != actual.rows() || m_expected.cols() != actual.cols())
                {
                    return false;
                }
                return actual.isApprox(m_expected, m_tolerance);
            }

            std::string describe() const override
            {
                std::ostringstream ss;
                ss << "is not approximately equal to:\n"
                   << m_expected
                   << "\nwithin tolerance " << m_tolerance;
                return ss.str();
            }

        private:
            const EigenType &m_expected;
            double m_tolerance;
        };

        template <typename EigenType>
        EigenApproxMatcher<EigenType> IsApprox(const EigenType &expected, double tolerance = 1e-9)
        {
            return EigenApproxMatcher<EigenType>(expected, tolerance);
        }
    }
}

#endif // GALILEO_TESTING_UTILS_MATCHERS_EIGEN_MATCHERS_HPP
