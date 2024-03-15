#pragma once

#include <pinocchio/autodiff/casadi.hpp>
#include <vector>
#include <string>

namespace galileo
{
    namespace opt
    {
        /**
         * @brief Helper class for storing polynomial information.
         *
         */
        class LagrangePolynomial
        {
        public:
            /**
             * @brief Construct a new Lagrange Polynomial object.
             *
             */
            LagrangePolynomial(){};

            /**
             * @brief Construct a new Lagrange Polynomial object. Compute and store the coefficients for a given degree and collocation scheme.
             * Use d_ = 0 for a piecewise constant polynomial.
             *
             * @param d_ Degree of the polynomial
             * @param scheme Collocation scheme: "radau" or "legendre"
             */
            LagrangePolynomial(int d_, const std::string &scheme = "radau");

            /**
             * @brief Perform symbolic Lagrange Interpolation, which, given a time from the Lagrange time scale, interpolates terms to find the value at time t.
             *
             * @param t Time to interpolate at
             * @param terms Terms at knot points to use for interpolation
             * @return Scalar Interpolated value
             */
            template <typename Scalar>
            Scalar barycentricInterpolation(double t, const std::vector<Scalar> terms) const;

            /**
             * @brief Perform symbolic Lagrange Interpolation, which, given a time from the Lagrange time scale, interpolates terms to find the value at time t.
             *
             * @param t Time to interpolate at
             * @param terms Terms at knot points to use for interpolation
             * @return Scalar Interpolated value
             */
            Eigen::VectorXd barycentricInterpolation(double t, const Eigen::MatrixXd &terms) const;

            /**
             * @brief If true, interprets the polynomial as a piecewise constant value
             * 
             */
            bool piecewise_constant = false;

            /**
             * @brief Degree of the polynomial.
             *
             */
            int d;

            /**
             * @brief The roots of the polynomial.
             *
             */
            std::vector<double> tau_root;

            /**
             * @brief Quadrature coefficients.
             *
             */
            std::vector<double> B;

            /**
             * @brief Collocation coeffficients.
             *
             */
            std::vector<std::vector<double>> C;

            /**
             * @brief Continuity coefficients.
             *
             */
            std::vector<double> D;

            /**
             * @brief Weights used for barycentric interpolation.
             * 
             */
            std::vector<double> barycentric_weights;
        };
    }
}