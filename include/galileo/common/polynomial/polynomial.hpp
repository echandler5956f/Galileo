#ifndef __galileo_common_polynomial_polynomial_hpp__
#define __galileo_common_polynomial_polynomial_hpp__

#include "galileo/common/fwd.hpp"
#include "galileo/common/polynomial/helpers.hpp"

namespace galileo
{
    namespace math
    {

        template <typename _NumScalar, int _N, int _Options>
        class JacobiPolynomialTpl
        {
        public:
            using NumScalar = _NumScalar;
            static constexpr int N = _N;
            static constexpr int Options = _Options;

            JacobiPolynomialTpl() : alpha_(0), beta_(0) {}

            JacobiPolynomialTpl(const NumScalar &alpha, const NumScalar &beta) : alpha_(alpha), beta_(beta)
            {
                compute_nodes_and_weights();
                compute_coefficients();
                compute_barycentric_weights();
            }

            template <typename InputMatrixType, typename OutputVectorType>
            void barycentricInterpolation(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &w, Eigen::MatrixBase<OutputVectorType> &u) const
            {
                typedef typename Eigen::internal::plain_row_type<OutputVectorType>::type RowVectorType;

                assert(w.cols() == N);
                assert(t >= 0. && t <= 1.);

                // Compute the interpolated value
                RowVectorType numerator = RowVectorType::Zero(w.rows());
                RowVectorType denominator = RowVectorType::Zero(w.rows());
                NumScalar interpolant;
                for (std::size_t i = 0; i < N; ++i)
                {
                    if (std::abs(t - nodes_[i]) < 1e-8)
                    {
                        u = w.col(i);
                        return;
                    }
                    interpolant = barycentric_weights_(i) / (t - nodes_(i));
                    numerator += interpolant * w.col(i);
                    denominator += RowVectorType::Constant(w.rows(), interpolant);
                }

                if ((denominator.array() == 0).any())
                {
                    throw std::runtime_error("Error: Division by zero in BarycentricInterpolation");
                }
                u = numerator.array() / denominator.array();
            }

            template <typename InputMatrixType, typename OutputMatrixType>
            void barycentricInterpolationDiff(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &w, Eigen::MatrixBase<OutputMatrixType> &du_dw) const
            {
                assert(w.cols() == N);
                assert(t >= NumScalar(0.) && t <= NumScalar(1.));

                assert(du_dw.rows() == w.rows());
                assert(du_dw.cols() == w.rows() * N);

                // If t is very close to one of the nodes, the interpolation directly returns w.col(i).
                // In that case, the sensitivity with respect to that column is the identity,
                // and with respect to all other columns is zero.
                for (std::size_t i = 0; i < N; ++i)
                {
                    if (std::abs(t - nodes_[i]) < 1e-8)
                    {
                        du_dw.setZero();
                        du_dw.block(0, i * w.rows(), w.rows(), w.rows()) = Eigen::Matrix<NumScalar, w.rows(), w.rows, Options>::Identity();
                        return;
                    }
                }

                // Compute the barycentric coefficients c_i and their sum.
                Eigen::Matrix<NumScalar, N, 1, Options> c;
                c.setZero();
                NumScalar sum_c = 0.0;
                for (std::size_t i = 0; i < N; ++i)
                {
                    c[i] = barycentric_weights_[i] / (t - nodes_[i]);
                    sum_c += c[i];
                }

                // Check for division by zero
                if (std::abs(sum_c) < 1e-12)
                {
                    throw std::runtime_error("Error: Division by zero in barycentricInterpolationDiff");
                }

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

                Eigen::Matrix<NumScalar, N, 1, Options> bj;
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
                imtqlx(nodes_, bj, weights_);

                // Map nodes_ to [0, 1]
                nodes_ = (nodes_ + 1.0) / 2.0;

                for (int i = 0; i < N; i++)
                {
                    weights_(i) = weights_(i) * weights_(i);
                }
            }

            void compute_coefficients()
            {
                // Diagonal matrix of nodes_
                Eigen::Matrix<NumScalar, N, N, Options> nodes_diag = nodes_.asDiagonal();

                // R is a matrix defined by diag(1. / (1 : N))
                Eigen::Matrix<NumScalar, N, 1, Options> R;
                for (int j = 0; j < N; j++)
                {
                    R(j, 0) = 1. / (j + 1);
                }
                Eigen::Matrix<NumScalar, N, N, Options> R_diag = R.asDiagonal();
                Eigen::Matrix<NumScalar, N, N, Options> Vandermonde = Eigen::Matrix<NumScalar, N, N, Options>::Ones();
                // Vandermonde matrix =
                //[[1, nodes_1, nodes_1 ^ 2, ..., nodes_1 ^(N - 1)],
                //[1, nodes_2, nodes_2 ^ 2, ..., nodes_2 ^(N - 1)],
                //...
                //[1, nodes_N, nodes_N ^ 2, ..., nodes_N ^(N - 1)]]
                for (int j = 0; j < N; j++)
                {
                    for (int r = 0; r < N; r++)
                    {
                        Vandermonde(j, r) = std::pow(nodes_[j + 1], r);
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

            Eigen::Matrix<NumScalar, N, N, Options> coeffs_;
            Eigen::Matrix<NumScalar, N, 1, Options> weights_;
            Eigen::Matrix<NumScalar, N, 1, Options> nodes_;

            Eigen::Matrix<NumScalar, N, 1, Options> barycentric_weights_;

            NumScalar alpha_;
            NumScalar beta_;

        }; // class JacobiPolynomialTpl

    } // namespace math

} // namespace galileo

#endif // __galileo_common_polynomial_polynomial_hpp__