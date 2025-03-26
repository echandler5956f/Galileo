#ifndef __galileo_math_poly_helpers_hpp__
#define __galileo_math_poly_helpers_hpp__

#include "galileo/math/fwd.hpp"


namespace galileo
{
    namespace math
    {

        // Inline specialized sign function.
        template <typename Scalar>
        inline Scalar r8_sign(Scalar x) noexcept
        {
            return (x < Scalar(0)) ? Scalar(-1) : Scalar(1);
        }

        // IMTQLX: Diagonalizes a symmetric tridiagonal matrix.
        template <typename VectorType1, typename VectorType2, typename VectorType3>
        inline void imtqlx(Eigen::MatrixBase<VectorType1> &d,
                           Eigen::MatrixBase<VectorType2> &e,
                           Eigen::MatrixBase<VectorType3> &z)
        {
            using Scalar = typename VectorType1::Scalar;
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

    } // namespace math

} // namespace galileo

#endif // __galileo_math_poly_helpers_hpp__