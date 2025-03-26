#ifndef __galileo_math_polynomial_hpp__
#define __galileo_math_polynomial_hpp__

#include "galileo/math/fwd.hpp"
#include "galileo/math/poly-helpers.hpp"

namespace galileo
{
    namespace math
    {

        template <typename _NumScalar, int _N, int _Options>
        class JacobiPolynomial
        {
        public:
            using NumScalar = _NumScalar;
            static constexpr int N = _N;
            static constexpr int Options = _Options;

            JacobiPolynomial(const NumScalar &alpha, const NumScalar &beta) : alpha_(alpha), beta_(beta)
            {
                compute_nodes_and_weights();
                compute_coefficients();
                compute_barycentric_weights();
            }

            template <typename InputMatrixType, typename OutputVectorType>
            inline void BarycentricInterpolation(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &terms, Eigen::MatrixBase<OutputVectorType> &output) const
            {
                typedef typename Eigen::internal::plain_row_type<OutputVectorType>::type RowVectorType;

                assert(terms.cols() == N);
                assert(t >= -1e-8 && t <= 1. + 1e-8);

                // Compute the interpolated value
                RowVectorType numerator = RowVectorType::Zero(terms.rows());
                RowVectorType denominator = RowVectorType::Zero(terms.rows());
                NumScalar interpolant;
                for (std::size_t i = 0; i < N; ++i)
                {
                    if (std::abs(t - nodes_[i]) < 1e-6)
                    {
                        output = terms.col(i);
                        return;
                    }
                    interpolant = barycentric_weights_(i) / (t - nodes_(i));
                    numerator += interpolant * terms.col(i);
                    denominator += RowVectorType::Constant(terms.rows(), interpolant);
                }

                if ((denominator.array() == 0).any())
                {
                    throw std::runtime_error("Error: Division by zero in BarycentricInterpolation");
                }
                output = numerator.array() / denominator.array();
            }

        protected:
            void compute_nodes_and_weights()
            {
                NumScalar ab = alpha + beta;
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
                    NumScalar i_r8 = Scalar(i + 1);
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
                Eigen::Matrix<NumScalar, N, N, Options> R = R.asDiagonal();
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
                coeffs_ = nodes_diag * Vandermonde * R * Vandermonde.inverse();
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

        }; // class JacobiPolynomial

    } // namespace math

} // namespace galileo

#endif // __galileo_math_polynomial_hpp__