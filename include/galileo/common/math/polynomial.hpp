#ifndef __galileo_common_math_polynomial_hpp__
#define __galileo_common_math_polynomial_hpp__

#include "galileo/common/fwd.hpp"
#include <cassert>
#include <limits>

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
        // d: diagonal elements
        // e: off-diagonal elements
        // z: on input, a vector. On output, the value of Q^T * z where Q is the orthogonal matrix that diagonalizes the matrix.
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

                    assert(j < itn && "IMTQLX - Fatal error!");
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

        template <typename _NumScalar, int _N, int _Options>
        class JacobiPolynomialTpl
        {
        public:
            using NumScalar = _NumScalar;
            static constexpr int N = _N;
            static constexpr int Options = _Options;

            JacobiPolynomialTpl() : alpha_(0.0), beta_(0.0) {}

            JacobiPolynomialTpl(const NumScalar &alpha, const NumScalar &beta) : alpha_(alpha), beta_(beta)
            {
                compute_nodes_and_weights();
                compute_coefficients();
                compute_barycentric_weights();
            }

            const Eigen::GMatrix<NumScalar, N, N, Options> &get_coeffs() const
            {
                return coeffs_;
            }

            const Eigen::GMatrix<NumScalar, N, 1, Options> &get_weights() const
            {
                return weights_;
            }

            const Eigen::GMatrix<NumScalar, N, 1, Options> &get_nodes() const
            {
                return nodes_;
            }

            const Eigen::GMatrix<NumScalar, N, N, Options> &get_jacobi_matrix() const
            {
                return jacobi_matrix_;
            }

            // Interpolate a polynomial at t with values w at the nodes, yielding the vector u(t)
            template <typename InputMatrixType, typename OutputVectorType>
            void barycentricInterpolation(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &w, Eigen::MatrixBase<OutputVectorType> &u) const
            {
                assert(w.cols() == N);
                assert(t >= NumScalar(0.0) && t <= NumScalar(1.0));

                for (std::size_t i = 0; i < N; ++i)
                {
                    if (std::abs(t - nodes_(i)) < std::numeric_limits<NumScalar>::epsilon())
                    {
                        u = w.col(i);
                        return;
                    }
                }

                Eigen::GMatrix<NumScalar, N, 1, Options> c;
                NumScalar sum_c = 0.0;
                for (std::size_t i = 0; i < N; ++i)
                {
                    c(i) = barycentric_weights_(i) / (t - nodes_(i));
                    sum_c += c(i);
                }

                assert(std::abs(sum_c) > std::numeric_limits<NumScalar>::epsilon() && "Error: Division by zero in BarycentricInterpolation");

                u.setZero();
                for (std::size_t i = 0; i < N; ++i)
                {
                    u += c(i) * w.col(i);
                }
                u /= sum_c;
            }

            template <typename InputMatrixType, typename OutputMatrixType>
            void barycentricInterpolationDiff(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &w, Eigen::MatrixBase<OutputMatrixType> &du_dw) const
            {
                assert(w.cols() == N);
                assert(t >= NumScalar(0.0) && t <= NumScalar(1.0));

                assert(du_dw.rows() == w.rows());
                assert(du_dw.cols() == w.rows() * N);

                // If t is very close to one of the nodes, the interpolation directly returns w.col(i).
                // In that case, the sensitivity with respect to that column is the identity,
                // and with respect to all other columns is zero.
                for (std::size_t i = 0; i < N; ++i)
                {
                    if (std::abs(t - nodes_[i]) < std::numeric_limits<NumScalar>::epsilon())
                    {
                        du_dw.setZero();
                        du_dw.block(0, i * w.rows(), w.rows(), w.rows()).setIdentity();
                        return;
                    }
                }

                // Compute the barycentric coefficients c_i and their sum.
                Eigen::GMatrix<NumScalar, N, 1, Options> c;
                c.setZero();
                NumScalar sum_c = 0.0;
                for (std::size_t i = 0; i < N; ++i)
                {
                    c[i] = barycentric_weights_[i] / (t - nodes_[i]);
                    sum_c += c[i];
                }

                // Check for division by zero
                assert(std::abs(sum_c) > std::numeric_limits<NumScalar>::epsilon() && "Error: Division by zero in barycentricInterpolationDiff");

                // The interpolated value is u = (sum_i c_i * w.col(i)) / sum_c.
                // Thus, for each element k of u:
                //     u(k) = (sum_i c_i * w(k,i)) / sum_c.
                // Therefore, the partial derivative with respect to w(k,j) is c_j/sum_c (if the row index matches).
                // We pack these derivatives into a vector of matrices, where each matrix is (w.rows() x w.rows())
                // representing the derivative with respect to one column of w.
                for (std::size_t j = 0; j < N; ++j)
                {
                    // For each column j, the sensitivity matrix is diagonal with constant c[j] / sum_c.
                    du_dw.block(0, j * w.rows(), w.rows(), w.rows()).diagonal().setConstant(c[j] / sum_c);
                }
            }

        protected:
            void compute_nodes_and_weights()
            {
                NumScalar ab = alpha_ + beta_;
                NumScalar abi = 2.0 + ab;

                // Define the zero-th moment.
                NumScalar zemu = std::pow(2.0, ab + 1.0) * std::tgamma(alpha_ + 1.0) * std::tgamma(beta_ + 1.0) / std::tgamma(abi);

                // Define the Jacobi matrix.
                nodes_(0) = (beta_ - alpha_) / abi;
                for (int i = 1; i < N; i++)
                {
                    nodes_(i) = 0.0;
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
                    nodes_(i) = a2b2 / ((abi - 2.0) * abi);
                    abi = abi * abi;
                    bj(i) = 4.0 * i_r8 * (i_r8 + alpha_) * (i_r8 + beta_) * (i_r8 + ab) / ((abi - 1.0) * abi);
                }

                for (int i = 0; i < N; i++)
                {
                    bj(i) = std::sqrt(bj(i));
                }

                weights_(0) = std::sqrt(zemu);
                for (int i = 1; i < N; i++)
                {
                    weights_(i) = 0.0;
                }

                // Diagonalize the Jacobi matrix.
                jacobi_matrix_.setZero();
                jacobi_matrix_.diagonal() = nodes_;
                jacobi_matrix_.diagonal(-1) = bj.segment(1, N - 1);
                jacobi_matrix_.diagonal(1) = bj.segment(1, N - 1);

                imtqlx(nodes_, bj, weights_);

                // Map nodes_ to [0, 1]
                nodes_ = (nodes_.array() + 1.0) / 2.0;

                for (int i = 0; i < N; i++)
                {
                    weights_(i) = weights_(i) * weights_(i);
                }
            }

            void compute_coefficients()
            {
                // Diagonal matrix of nodes_
                Eigen::GMatrix<NumScalar, N, N, Options> nodes_diag = nodes_.asDiagonal();

                // R is a matrix defined by diag(1. / (1 : N))
                Eigen::GMatrix<NumScalar, N, 1, Options> R;
                for (int j = 0; j < N; j++)
                {
                    R(j, 0) = 1.0 / (j + 1);
                }
                Eigen::GMatrix<NumScalar, N, N, Options> R_diag = R.asDiagonal();
                Eigen::GMatrix<NumScalar, N, N, Options> Vandermonde = Eigen::GMatrix<NumScalar, N, N, Options>::Ones();
                // Vandermonde matrix =
                //[[1, nodes_1, nodes_1 ^ 2, ..., nodes_1 ^(N - 1)],
                //[1, nodes_2, nodes_2 ^ 2, ..., nodes_2 ^(N - 1)],
                //...
                //[1, nodes_N, nodes_N ^ 2, ..., nodes_N ^(N - 1)]]
                for (int j = 0; j < N; j++)
                {
                    for (int r = 0; r < N; r++)
                    {
                        Vandermonde(j, r) = std::pow(nodes_(j), r);
                    }
                }

                // Only working for nodes \in[0, 1]
                coeffs_ = nodes_diag * Vandermonde * R_diag * Vandermonde.inverse();
            }

            void compute_barycentric_weights()
            {
                /*Barycentric weights*/
                for (int j = 0; j < N; ++j)
                    barycentric_weights_(j) = 1.0;

                /*For all collocation points*/
                for (int j = 0; j < N; ++j)
                {
                    for (int r = 0; r < N; ++r)
                    {
                        if (r != j)
                        {
                            barycentric_weights_(j) *= (nodes_(j) - nodes_(r));
                        }
                    }
                    barycentric_weights_(j) = 1.0 / barycentric_weights_(j);
                }
            }

            Eigen::GMatrix<NumScalar, N, N, Options> coeffs_;
            Eigen::GMatrix<NumScalar, N, 1, Options> weights_;
            Eigen::GMatrix<NumScalar, N, 1, Options> nodes_;
            Eigen::GMatrix<NumScalar, N, N, Options> jacobi_matrix_;

            Eigen::GMatrix<NumScalar, N, 1, Options> barycentric_weights_;

            NumScalar alpha_;
            NumScalar beta_;

        }; // class JacobiPolynomialTpl

    } // namespace math

} // namespace galileo

#endif // __galileo_common_math_polynomial_hpp__
