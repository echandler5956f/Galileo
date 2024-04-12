#include "galileo/opt/LagrangePolynomial.h"

namespace galileo
{
    namespace opt
    {
        LagrangePolynomial::LagrangePolynomial(int d_, const std::string &scheme)
        {
            if (d_ == 0)
            {
                piecewise_constant = true;
                d = 1;
            }
            else
                d = d_;
            /*Choose collocation points*/
            /*IMPORTANT NOTE for radau scheme: The right endpoint is NOT included in the collocation points.*/
            tau_root = casadi::collocation_points(d, scheme);
            tau_root.insert(tau_root.begin(), 0);

            /*Coefficients of the quadrature function*/
            B.resize(d + 1);

            /*Coefficients of the collocation equation*/
            C.resize(d + 1);
            for (int j = 0; j < d + 1; ++j)
                C[j].resize(d + 1);

            /*Coefficients of the continuity equation*/
            D.resize(d + 1);

            /*Barycentric weights*/
            barycentric_weights.resize(d + 1);
            for (int j = 0; j < d + 1; ++j)
                barycentric_weights[j] = 1.0;

            /*For all collocation points*/
            for (int j = 0; j < d + 1; ++j)
            {
                /*Construct Lagrange polynomials to get the polynomial basis at the collocation point*/
                casadi::Polynomial p = 1;
                for (int r = 0; r < d + 1; ++r)
                {
                    if (r != j)
                    {
                        p *= casadi::Polynomial(-tau_root[r], 1) / (tau_root[j] - tau_root[r]);
                        barycentric_weights[j] *= (tau_root[j] - tau_root[r]);
                    }
                }
                barycentric_weights[j] = 1.0 / barycentric_weights[j];
                /*Evaluate the polynomial at the final time to get the coefficients of the continuity equation*/
                D[j] = p(1.0);

                /*Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation*/
                casadi::Polynomial dp = p.derivative();
                for (int r = 0; r < d + 1; ++r)
                {
                    C[j][r] = dp(tau_root[r]);
                }
                casadi::Polynomial pint = p.anti_derivative();
                B[j] = pint(1.0);
            }
        }

        template <typename Scalar>
        Scalar LagrangePolynomial::barycentricInterpolation(double t, const std::vector<Scalar> terms) const
        {
            assert(terms.size() == tau_root.size());
            assert(t >= -1e-8 && t <= 1. + 1e-8);
            if (piecewise_constant)
            {
                return terms[0];
            }

            // Compute the interpolated value
            Scalar numerator = Scalar::zeros(terms[0].size1(), 1);
            Scalar denominator = Scalar::zeros(terms[0].size1(), 1);
            Scalar interpolant;
            for (std::size_t i = 0; i < tau_root.size(); ++i)
            {
                if (std::abs(t - tau_root[i]) < 1e-6)
                {
                    return terms[i];
                }
                interpolant = Scalar(barycentric_weights[i]) / (t - tau_root[i]);
                numerator += interpolant * terms[i];
                denominator += interpolant;
            }

            return numerator / denominator;
        }

        template casadi::DM LagrangePolynomial::barycentricInterpolation<casadi::DM>(double t, const std::vector<casadi::DM> terms) const;
        template casadi::SX LagrangePolynomial::barycentricInterpolation<casadi::SX>(double t, const std::vector<casadi::SX> terms) const;
        template casadi::MX LagrangePolynomial::barycentricInterpolation<casadi::MX>(double t, const std::vector<casadi::MX> terms) const;

        Eigen::VectorXd LagrangePolynomial::barycentricInterpolation(double t, const Eigen::MatrixXd &terms) const
        {
            assert(terms.cols() == tau_root.size());
            assert(t >= -1e-8 && t <= 1. + 1e-8);
            if (piecewise_constant)
            {
                return terms.col(0);
            }

            // Compute the interpolated value
            Eigen::VectorXd numerator = Eigen::VectorXd::Zero(terms.rows());
            Eigen::VectorXd denominator = Eigen::VectorXd::Zero(terms.rows());
            double interpolant;
            for (std::size_t i = 0; i < tau_root.size(); ++i)
            {
                if (std::abs(t - tau_root[i]) < 1e-6)
                {
                    return terms.col(i);
                }
                interpolant = barycentric_weights[i] / (t - tau_root[i]);
                numerator += interpolant * terms.col(i);
                denominator += Eigen::VectorXd::Constant(terms.rows(), interpolant);
            }

            if ((denominator.array() == 0).any())
            {
                throw std::runtime_error("Error: Division by zero in barycentricInterpolation");
            }
            return numerator.array() / denominator.array();
        }
    }
}