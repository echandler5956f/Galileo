#ifndef GALILEO_TESTING_UTILS_MATCHERS_FINITE_DIFFERENCE_HPP
#define GALILEO_TESTING_UTILS_MATCHERS_FINITE_DIFFERENCE_HPP

#include <functional>
#include <string>
#include <Eigen/Core>
#include "galileo/core/states/state-base.hpp"
#include "tests/utils/matchers/eigen_matchers.hpp"

namespace galileo
{
    namespace testing
    {
        /**
         * @brief Checks a Jacobian matrix by comparing it against a numerical one computed with manifold-aware finite differences.
         *
         * This function is the core of derivative checking in Galileo. It perturbs the input vector `x` on its manifold
         * using the `state.integrate` method and computes the change in the output of `func` using `state.diff`.
         * This correctly handles the geometry of the state space.
         *
         * @tparam StateType The state model type, derived from StateBase.
         * @tparam Func The function whose Jacobian is to be tested. It should take an Eigen vector and return an Eigen vector.
         * @tparam AnalyticalJacobian The type of the analytically computed Jacobian matrix.
         * @param state The state object that provides the manifold operations.
         * @param J_analytical The analytically computed Jacobian to be verified.
         * @param func The function to differentiate.
         * @param x The point at which to compute the Jacobian.
         * @param test_name A descriptive name for the test section.
         * @param tolerance The tolerance for the IsApprox check.
         * @param eps The step size for the finite difference approximation.
         */
        template <typename StateType, typename Func, typename AnalyticalJacobian>
        void check_jacobian_manifold(const StateType &state,
                                     const Eigen::MatrixBase<AnalyticalJacobian> &J_analytical,
                                     const Func &func,
                                     const typename StateType::VectorNx_t &x,
                                     const std::string &test_name,
                                     double tolerance = 1e-5, double eps = 1e-7)
        {
            using VectorNx_t = typename StateType::VectorNx_t;
            using VectorNdx_t = typename StateType::VectorNdx_t;

            auto y = func(x);
            AnalyticalJacobian J_numerical(y.size(), state.get_ndx());
            VectorNdx_t dx_perturb = VectorNdx_t::Zero(state.get_ndx());
            VectorNx_t x_perturbed;

            for (int i = 0; i < state.get_ndx(); ++i)
            {
                dx_perturb(i) = eps;
                state.integrate(x, dx_perturb, x_perturbed);
                auto y_plus = func(x_perturbed);

                dx_perturb(i) = -eps;
                state.integrate(x, dx_perturb, x_perturbed);
                auto y_minus = func(x_perturbed);

                J_numerical.col(i) = (y_plus - y_minus) / (2.0 * eps);
                dx_perturb(i) = 0.0;
            }

            CAPTURE(test_name);
            REQUIRE_THAT(J_analytical, IsApprox(J_numerical, tolerance));
        }

        /**
         * @brief Overload for check_jacobian_manifold for when the variable being differentiated is in a Euclidean space (e.g., a tangent vector), but the function still depends on a manifold state.
         */
        template <typename StateType, typename Func, typename InputType, typename AnalyticalJacobian>
        void check_jacobian_manifold(const StateType &state,
                                     const Eigen::MatrixBase<AnalyticalJacobian> &J_analytical,
                                     const Func &func,
                                     const Eigen::MatrixBase<InputType> &x_euclidean,
                                     const std::string &test_name,
                                     double tolerance = 1e-5, double eps = 1e-7)
        {
            auto y = func(x_euclidean);
            AnalyticalJacobian J_numerical(y.size(), x_euclidean.size());
            InputType x_perturbed = x_euclidean;

            for (int i = 0; i < x_euclidean.size(); ++i)
            {
                double h = eps;
                x_perturbed(i) += h;
                auto y_plus = func(x_perturbed);

                x_perturbed(i) -= 2 * h;
                auto y_minus = func(x_perturbed);

                J_numerical.col(i) = (y_plus - y_minus) / (2.0 * h);
                x_perturbed(i) = x_euclidean(i);
            }

            CAPTURE(test_name);
            REQUIRE_THAT(J_analytical, IsApprox(J_numerical, tolerance));
        }

        /**
         * @brief Checks a Jacobian matrix for a simple Euclidean function.
         */
        template <typename Func, typename InputType, typename AnalyticalJacobian>
        void check_jacobian_euclidean(const Eigen::MatrixBase<AnalyticalJacobian> &J_analytical,
                                      const Func &func,
                                      const Eigen::MatrixBase<InputType> &x,
                                      const std::string &test_name,
                                      double tolerance = 1e-5, double eps = 1e-7)
        {
            auto y = func(x);
            AnalyticalJacobian J_numerical(y.size(), x.size());
            InputType x_perturbed = x;

            for (int i = 0; i < x.size(); ++i)
            {
                double h = eps;
                x_perturbed(i) += h;
                auto y_plus = func(x_perturbed);

                x_perturbed(i) -= 2 * h;
                auto y_minus = func(x_perturbed);

                J_numerical.col(i) = (y_plus - y_minus) / (2.0 * h);
                x_perturbed(i) = x(i);
            }

            CAPTURE(test_name);
            REQUIRE_THAT(J_analytical, IsApprox(J_numerical, tolerance));
        }
    }
}

#endif // GALILEO_TESTING_UTILS_MATCHERS_FINITE_DIFFERENCE_HPP
