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
                    }
                }
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
            if (piecewise_constant)
            {
                return terms[0];
            }
            else
            {
                /*Sometimes we get negatives here, but I don't think we should since tau_root is [0,1]? Yet, the interpolation seems to be correct..*/
                // std::cout << "Interpolating at time " << t << std::endl;
                std::vector<Scalar> w;
                for (std::size_t j = 0; j < tau_root.size(); ++j)
                {
                    w.push_back(1.0);
                }
                // Compute the barycentric weights
                for (int j = 0; j < d; ++j)
                {
                    for (std::size_t r = 0; r < tau_root.size(); ++r)
                    {
                        if (r != std::size_t(j))
                        {
                            w[j] *= (tau_root[j] - tau_root[r]);
                        }
                    }
                    w[j] = 1.0 / w[j];
                }
                // Compute the interpolated value
                Scalar numerator = 0.0, denominator = 0.0;
                for (std::size_t i = 0; i < tau_root.size(); ++i)
                {
                    if (std::abs(t - tau_root[i]) < 1e-6)
                    {
                        return terms[i];
                    }
                    Scalar term = w[i] / (t - tau_root[i]);
                    numerator += term * terms[i];
                    denominator += term;
                }

                return numerator / denominator;
            }
        }

        template double LagrangePolynomial::barycentricInterpolation<double>(double t, const std::vector<double> terms) const;
        template casadi::DM LagrangePolynomial::barycentricInterpolation<casadi::DM>(double t, const std::vector<casadi::DM> terms) const;
        template casadi::SX LagrangePolynomial::barycentricInterpolation<casadi::SX>(double t, const std::vector<casadi::SX> terms) const;
        template casadi::MX LagrangePolynomial::barycentricInterpolation<casadi::MX>(double t, const std::vector<casadi::MX> terms) const;
    }
}