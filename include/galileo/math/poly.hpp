#ifndef __galileo_math_poly_hpp__
#define __galileo_math_poly_hpp__

#include "galileo/math/fwd.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <tgmath.h>

namespace galileo
{
    namespace math
    {

        //----------------------------------------------------------------------
        // Inline specialized sign function.
        //----------------------------------------------------------------------
        template <typename Scalar>
        inline Scalar r8_sign(Scalar x) noexcept
        {
            return (x < Scalar(0)) ? Scalar(-1) : Scalar(1);
        }

        template <typename VectorType>
        inline void imtqlx(Eigen::MatrixBase<VectorType> &d,
                           Eigen::MatrixBase<VectorType> &e,
                           Eigen::MatrixBase<VectorType> &z)
        {
            using Scalar = typename VectorType::Scalar;
            const std::size_t n = d.size();

            Scalar b;
            Scalar c;
            Scalar f;
            Scalar g;
            int i;
            int ii;
            int itn = 30;
            int j;
            int k;
            int l;
            int m;
            int mml;
            Scalar p;
            Scalar prec;
            Scalar r;
            Scalar s;

            prec = std::numeric_limits<Scalar>::epsilon();

            if (n == 1)
            {
                return;
            }

            e(n - 1) = 0.0;

            for (l = 1; l <= n; l++)
            {
                j = 0;
                for (;;)
                {
                    for (m = l; m <= n; m++)
                    {
                        if (m == n)
                        {
                            break;
                        }
                        if (std::fabs(e(m - 1)) <= prec * (std::fabs(d(m - 1)) + std::fabs(d(m))))
                        {
                            break;
                        }
                    }
                    p = d(l - 1);
                    if (m == l)
                    {
                        break;
                    }

                    assert(itn <= j && "IMTQLX - Fatal error!");
                    j = j + 1;
                    g = (d(l) - p) / (2.0 * e(l - 1));
                    r = std::sqrt(g * g + 1.0);
                    g = d(m - 1) - p + e(l - 1) / (g + std::fabs(r) * r8_sign(g));
                    s = 1.0;
                    c = 1.0;
                    p = 0.0;
                    mml = m - l;

                    for (ii = 1; ii <= mml; ii++)
                    {
                        i = m - ii;
                        f = s * e(i - 1);
                        b = c * e(i - 1);

                        if (std::fabs(g) <= std::fabs(f))
                        {
                            c = g / f;
                            r = std::sqrt(c * c + 1.0);
                            e(i) = f * r;
                            s = 1.0 / r;
                            c = c * s;
                        }
                        else
                        {
                            s = f / g;
                            r = std::sqrt(s * s + 1.0);
                            e(i) = g * r;
                            c = 1.0 / r;
                            s = s * c;
                        }
                        g = d(i) - p;
                        r = (d(i - 1) - g) * s + 2.0 * c * b;
                        p = s * r;
                        d(i) = g + p;
                        g = c * r - b;
                        f = z(i);
                        z(i) = s * z(i - 1) + c * f;
                        z(i - 1) = c * z(i - 1) - s * f;
                    }
                    d(l - 1) = d(l - 1) - p;
                    e(l - 1) = g;
                    e(m - 1) = 0.0;
                }
            }

            //  Sorting.
            for (ii = 2; ii <= m; ii++)
            {
                i = ii - 1;
                k = i;
                p = d(i - 1);

                for (j = ii; j <= n; j++)
                {
                    if (d(j - 1) < p)
                    {
                        k = j;
                        p = d(j - 1);
                    }
                }
                if (k != i)
                {
                    d(k - 1) = d(i - 1);
                    d(i - 1) = p;
                    p = z(i - 1);
                    z(i - 1) = z(k - 1);
                    z(k - 1) = p;
                }
            }
        }

        //==============================================================================
        // DYNAMIC-SIZED VERSION
        //==============================================================================

        template <typename Scalar>
        inline Eigen::Matrix<Scalar, -1, 1>
        j_polynomial_zeros(std::size_t n, Scalar alpha, Scalar beta)
        {
            using VectorType = Eigen::Matrix<Scalar, -1, 1>;

            if (n == 0)
            {
                throw std::invalid_argument("Order n must be positive.");
            }
            if (alpha <= Scalar(-1) || beta <= Scalar(-1))
            {
                throw std::invalid_argument("Parameters alpha and beta must be > -1.");
            }

            const Scalar ab = alpha + beta;
            const Scalar abi = Scalar(2) + ab;

            // Zero-th moment
            const Scalar zemu = std::pow(Scalar(2), ab + Scalar(1)) * std::tgamma(alpha + Scalar(1)) * std::tgamma(beta + Scalar(1)) / std::tgamma(abi);

            // Allocate the Jacobi matrix
            VectorType d(n), e(n);
            d.setZero();
            e.setZero();

            // First diagonal element, first subdiagonal element
            d(0) = (beta - alpha) / abi;
            e(0) = Scalar(4) * (Scalar(1) + alpha) * (Scalar(1) + beta) / ((abi + Scalar(1)) * abi * abi);

            const Scalar a2b2 = beta * beta - alpha * alpha;

            // Fill diagonal and subdiagonal
            for (std::size_t i = 1; i < n; ++i)
            {
                const Scalar i_r8 = Scalar(i + 1);
                const Scalar current_abi = Scalar(2) * i_r8 + ab;
                d(i) = a2b2 / ((current_abi - Scalar(2)) * current_abi);
                const Scalar next_abi = current_abi * current_abi;
                e(i) = Scalar(4) * i_r8 * (i_r8 + alpha) * (i_r8 + beta) * (i_r8 + ab) / ((next_abi - Scalar(1)) * next_abi);

                // Take sqrt in-place
                e(i) = std::sqrt(e(i));
            }
            // The first subdiagonal also needs sqrt
            e(0) = std::sqrt(e(0));

            // For IMTQLX, we also need z (the initial "vector" to be transformed).
            VectorType z(n);
            z.setZero();
            z(0) = std::sqrt(zemu);

            // Diagonalize the Jacobi matrix
            imtqlx(d, e, z);

            // The zeros are the eigenvalues stored in d
            return d;
        }

        //==============================================================================
        // FIXED-SIZE VERSION
        //   Allows for compiler optimizations (loop unrolling, etc.) when N is known
        //   at compile time.
        //==============================================================================

        template <typename Scalar, std::size_t N>
        inline Eigen::Matrix<Scalar, N, 1>
        j_polynomial_zeros(Scalar alpha, Scalar beta)
        {
            using VectorType = Eigen::Matrix<Scalar, N, 1>;
            constexpr std::size_t n = static_cast<std::size_t>(N);

            if (alpha <= Scalar(-1) || beta <= Scalar(-1))
            {
                throw std::invalid_argument("Parameters alpha and beta must be > -1.");
            }

            const Scalar ab = alpha + beta;
            const Scalar abi = Scalar(2) + ab;

            // Zero-th moment
            const Scalar zemu = std::pow(Scalar(2), ab + Scalar(1)) * std::tgamma(alpha + Scalar(1)) * std::tgamma(beta + Scalar(1)) / std::tgamma(abi);

            // Fixed-size allocations
            VectorType d, e, z;
            d.setZero();
            e.setZero();
            z.setZero();

            // First diagonal element, first subdiagonal
            d(0) = (beta - alpha) / abi;
            if constexpr (N > 1)
            {
                e(0) = Scalar(4) * (Scalar(1) + alpha) * (Scalar(1) + beta) / ((abi + Scalar(1)) * abi * abi);
                e(0) = std::sqrt(e(0));
            }

            const Scalar a2b2 = beta * beta - alpha * alpha;

            // Fill diagonal and subdiagonal
            for (std::size_t i = 1; i < n; ++i)
            {
                const Scalar i_r8 = Scalar(i + 1);
                const Scalar current_abi = Scalar(2) * i_r8 + ab;
                d(i) = a2b2 / ((current_abi - Scalar(2)) * current_abi);
                const Scalar next_abi = current_abi * current_abi;
                e(i) = Scalar(4) * i_r8 * (i_r8 + alpha) * (i_r8 + beta) * (i_r8 + ab) / ((next_abi - Scalar(1)) * next_abi);

                e(i) = std::sqrt(e(i));
            }

            // z for use in IMTQLX
            z(0) = std::sqrt(zemu);

            // Diagonalize
            imtqlx(d, e, z);

            // The zeros of the Jacobi polynomial are the eigenvalues (in d).
            return d;
        }

    } // namespace math

} // namespace galileo

#endif // __galileo_math_poly_hpp__