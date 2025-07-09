#ifndef __galileo_common_math_jacobi_roots_hpp__
#define __galileo_common_math_jacobi_roots_hpp__

#include "galileo/common/fwd.hpp"

#include <cassert>
#include <limits>

namespace galileo
{

    // IMTQLX: Diagonalizes a symmetric tridiagonal matrix.
    // d: diagonal elements
    // e: off-diagonal elements
    // z: on input, a vector. On output, the value of Q^T * z where Q is the orthogonal matrix that diagonalizes the matrix.
    template <typename VectorType1, typename VectorType2, typename VectorType3>
    inline void imtqlx(Eigen::MatrixBase<VectorType1> &d,
                       Eigen::MatrixBase<VectorType2> &e,
                       Eigen::MatrixBase<VectorType3> &z)
    {
        using Scalar = typename VectorType1::Scalar;
        const int n = d.size();

        Scalar b;
        Scalar c;
        Scalar f;
        Scalar g;
        Scalar sign_g;
        int i;
        int ii;
        const int itn = 30;
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

        for (l = 1; l <= (int)n; l++)
        {
            j = 0;
            for (;;)
            {
                for (m = l; m <= (int)n; m++)
                {
                    if (m == (int)n)
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

                GALILEO_ASSERT(j < itn, "IMTQLX: Fatal error!");
                j = j + 1;
                g = (d(l) - p) / (2.0 * e(l - 1));
                r = std::sqrt(g * g + 1.0);
                sign_g = (g < Scalar(0.)) ? Scalar(-1.) : Scalar(1.);
                g = d(m - 1) - p + e(l - 1) / (g + std::fabs(r) * sign_g);
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

            for (j = ii; j <= (int)n; j++)
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

    template <typename _NumScalar, int _N, int _Options = 0>
    class JacobiRootsTpl
    {
    public:
        using NumScalar = _NumScalar;
        static constexpr int N = _N;
        static constexpr int Options = _Options;

        JacobiRootsTpl() : alpha_(0.0), beta_(0.0) {}

        JacobiRootsTpl(const NumScalar &alpha, const NumScalar &beta) : alpha_(alpha), beta_(beta)
        {
        }

        void compute_roots()
        {
            NumScalar ab = alpha_ + beta_;
            NumScalar abi = 2.0 + ab;

            // Define the zero-th moment.
            NumScalar zemu = std::pow(2.0, ab + 1.0) * std::tgamma(alpha_ + 1.0) * std::tgamma(beta_ + 1.0) / std::tgamma(abi);

            // Define the Jacobi matrix.
            roots_(0) = (beta_ - alpha_) / abi;
            for (int i = 1; i < N; i++)
            {
                roots_(i) = 0.0;
            }

            Eigen::GMatrix<NumScalar, N, 1, Options> bj;
            bj.setZero();

            bj(0) = 4.0 * (1.0 + alpha_) * (1.0 + beta_) / ((abi + 1.0) * abi * abi);
            for (int i = 1; i < N; i++)
            {
                bj(i) = 0.0;
            }

            NumScalar a2b2 = beta_ * beta_ - alpha_ * alpha_;

            for (int i = 1; i < N; i++)
            {
                NumScalar i_r8 = NumScalar(i + 1);
                abi = 2.0 * i_r8 + ab;
                roots_(i) = a2b2 / ((abi - 2.0) * abi);
                abi = abi * abi;
                bj(i) = 4.0 * i_r8 * (i_r8 + alpha_) * (i_r8 + beta_) * (i_r8 + ab) / ((abi - 1.0) * abi);
            }

            for (int i = 0; i < N; i++)
            {
                bj(i) = std::sqrt(bj(i));
            }

            Eigen::GMatrix<NumScalar, N, 1, Options> gauss_jacobi_weights;
            gauss_jacobi_weights.setZero();

            gauss_jacobi_weights(0) = std::sqrt(zemu);
            for (int i = 1; i < N; i++)
            {
                gauss_jacobi_weights(i) = 0.0;
            }

            // Diagonalize the Jacobi matrix.
            jacobi_matrix_.setZero();
            jacobi_matrix_.diagonal() = roots_;
            jacobi_matrix_.diagonal(-1) = bj.segment(1, N - 1);
            jacobi_matrix_.diagonal(1) = bj.segment(1, N - 1);

            imtqlx(roots_, bj, gauss_jacobi_weights);

            // Map nodes_ to [0, 1]
            roots_ = (roots_.array() + 1.0) / 2.0;

            for (int i = 0; i < N; i++)
            {
                gauss_jacobi_weights(i) = gauss_jacobi_weights(i) * gauss_jacobi_weights(i);
            }
        }

        const Eigen::GMatrix<NumScalar, N, N, Options> &get_jacobi_matrix() const
        {
            return jacobi_matrix_;
        }

        const Eigen::GMatrix<NumScalar, N, 1, Options> &get_roots() const
        {
            return roots_;
        }

        NumScalar get_alpha() const
        {
            return alpha_;
        }

        NumScalar get_beta() const
        {
            return beta_;
        }

    protected:
        Eigen::GMatrix<NumScalar, N, N, Options> jacobi_matrix_;
        Eigen::GMatrix<NumScalar, N, 1, Options> roots_;

        NumScalar alpha_;
        NumScalar beta_;

    }; // class JacobiRootsTpl

} // namespace galileo

#endif // __galileo_common_math_jacobi_roots_hpp__
