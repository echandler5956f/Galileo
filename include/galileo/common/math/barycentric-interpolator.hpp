#ifndef __galileo_common_math_barycentric_interpolator_hpp__
#define __galileo_common_math_barycentric_interpolator_hpp__

#include "galileo/common/fwd.hpp"

#include <limits>

namespace galileo
{

    template <typename _NumScalar, int _N, int _Options = 0>
    class BarycentricInterpolatorTpl
    {
    public:
        using NumScalar = _NumScalar;
        static constexpr int N = _N;
        static constexpr int Options = _Options;

        using VectorN = Eigen::GMatrix<NumScalar, N, 1, Options>;
        using MatrixN = Eigen::GMatrix<NumScalar, N, N, Options>;

        template <typename InputVectorType>
        BarycentricInterpolatorTpl(const Eigen::MatrixBase<InputVectorType> &nodes)
        {
            // Set the dimension and validate consistency
            NDim_.set_value(nodes.size());

            // For fixed-size interpolators, check dimension consistency and copy element-wise
            if constexpr (DimensionTpl<N>::IsFixed)
            {
                GALILEO_ASSERT(nodes.size() == N,
                               "BarycentricInterpolator: Number of nodes must match template parameter N for "
                               "fixed-size interpolators.");
                if constexpr (InputVectorType::RowsAtCompileTime == Eigen::Dynamic ||
                              InputVectorType::ColsAtCompileTime == Eigen::Dynamic)
                {
                    nodes_ = nodes;
                }
                else
                {
                    // Element-wise copy to avoid compile-time dimension mismatch
                    for (int i = 0; i < N; ++i)
                    {
                        nodes_(i) = nodes(i);
                    }
                }
            }
            else
            {
                // For dynamic size, direct assignment is fine
                nodes_ = nodes;
            }

            compute_weights();
        }

        /**
         * @brief Compute interpolated values at time t
         *
         * @param t Interpolation parameter in [0, 1]
         * @param w Input matrix where each column represents the multidimensional state at a node
         *          Dimensions: (num_output_dimensions × num_nodes)
         * @param u Output vector containing the interpolated state
         *          Dimensions: (num_output_dimensions × 1)
         *
         * This method performs barycentric interpolation across multiple dimensions simultaneously.
         * Each column of w contains the values at one of the interpolation nodes, and the method
         * computes the weighted average to produce the interpolated values at time t.
         */
        template <typename InputMatrixType, typename OutputVectorType>
        void calc(const NumScalar &t,
                  const Eigen::MatrixBase<InputMatrixType> &w,
                  Eigen::MatrixBase<OutputVectorType> &u) const
        {
            using VectorType = typename Eigen::internal::plain_matrix_type<OutputVectorType>::type;

            GALILEO_ASSERT(
                w.cols() == NDim_.value(),
                "BarycentricInterpolator: Input matrix must have the same number of columns as the number of nodes.");
            GALILEO_ASSERT(t >= NumScalar(0.) && t <= NumScalar(1.),
                           "BarycentricInterpolator: Time must be between 0 and 1.");

            // Compute the interpolated value using proper barycentric formula
            VectorType numerator = VectorType::Zero(w.rows());
            NumScalar denominator = NumScalar(0.);
            NumScalar interpolant;
            for (int i = 0; i < NDim_.value(); ++i)
            {
                if (std::abs(t - nodes_(i)) < std::numeric_limits<NumScalar>::epsilon())
                {
                    u = w.col(i);
                    return;
                }
                interpolant = weights_(i) / (t - nodes_(i));
                numerator += interpolant * w.col(i);
                denominator += interpolant;
            }

            GALILEO_ASSERT(std::abs(denominator) > std::numeric_limits<NumScalar>::epsilon(),
                           "BarycentricInterpolator: Denominator is zero.");

            u = numerator / denominator;
        }

        /**
         * @brief Compute sensitivities of interpolated values with respect to input values
         *
         * @param t Interpolation parameter in [0, 1]
         * @param w Input matrix where each column represents the multidimensional state at a node
         *          Dimensions: (num_output_dimensions × num_nodes)
         * @param du_dw Output sensitivity matrix
         *          Dimensions: (num_output_dimensions × num_output_dimensions*num_nodes)
         *
         * The output matrix du_dw is organized as blocks, where each block of size
         * (num_output_dimensions × num_output_dimensions) represents the sensitivity of the
         * interpolated output with respect to one column of w. For barycentric interpolation,
         * these blocks are diagonal matrices.
         */
        template <typename InputMatrixType, typename OutputMatrixType>
        void calcDiff(const NumScalar &t,
                      const Eigen::MatrixBase<InputMatrixType> &w,
                      Eigen::MatrixBase<OutputMatrixType> &du_dw) const
        {
            GALILEO_ASSERT(
                w.cols() == NDim_.value(),
                "BarycentricInterpolator: Input matrix must have the same number of columns as the number of nodes.");
            GALILEO_ASSERT(t >= NumScalar(0.) && t <= NumScalar(1.),
                           "BarycentricInterpolator: Time must be between 0 and 1.");

            GALILEO_ASSERT(
                du_dw.rows() == w.rows(),
                "BarycentricInterpolator: Output matrix must have the same number of rows as the input matrix.");
            GALILEO_ASSERT(du_dw.cols() == w.rows() * NDim_.value(),
                           "BarycentricInterpolator: Output matrix must have (num_output_dims * num_nodes) columns.");

            // If t is very close to one of the nodes, the interpolation directly returns w.col(i).
            // In that case, the sensitivity with respect to that column is the identity,
            // and with respect to all other columns is zero.
            for (int i = 0; i < NDim_.value(); ++i)
            {
                if (std::abs(t - nodes_(i)) < std::numeric_limits<NumScalar>::epsilon())
                {
                    du_dw.setZero();
                    block(du_dw, 0, i * w.rows(), w.rows(), w.rows()) =
                        Eigen::Matrix<NumScalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(w.rows(), w.rows());
                    return;
                }
            }

            // Compute the barycentric coefficients c_i and their sum.
            VectorN c = VectorN::Zero(NDim_.value());
            NumScalar sum_c = 0.0;
            for (int i = 0; i < NDim_.value(); ++i)
            {
                c(i) = weights_(i) / (t - nodes_(i));
                sum_c += c(i);
            }

            // Check for division by zero
            GALILEO_ASSERT(std::abs(sum_c) > std::numeric_limits<NumScalar>::epsilon(),
                           "BarycentricInterpolator: Sum of coefficients is zero.");

            // The interpolated value is u = (sum_i c_i * w.col(i)) / sum_c.
            // Thus, for each element k of u:
            //     u(k) = (sum_i c_i * w(k,i)) / sum_c.
            // Therefore, the partial derivative with respect to w(k,j) is c_j/sum_c (if the row index matches).
            // We pack these derivatives into a vector of matrices, where each matrix is (w.rows() x w.rows())
            // representing the derivative with respect to one column of w.
            for (int j = 0; j < NDim_.value(); ++j)
            {
                // For each column j, the sensitivity matrix is diagonal with constant c[j] / sum_c.
                block(du_dw, 0, j * w.rows(), w.rows(), w.rows()).diagonal().setConstant(c(j) / sum_c);
            }
        }

        constexpr int get_n() const { return NDim_.value(); }

        const DimensionTpl<N> &NDim() const { return NDim_; }

    protected:
        void compute_weights()
        {
            weights_ = VectorN::Zero(NDim_.value());

            /*Barycentric weights*/
            for (int j = 0; j < NDim_.value(); ++j) weights_(j) = 1.0;

            /*For all nodes*/
            for (int j = 0; j < NDim_.value(); ++j)
            {
                for (int r = 0; r < NDim_.value(); ++r)
                {
                    if (r != j)
                    {
                        weights_(j) *= (nodes_(j) - nodes_(r));
                    }
                }
                weights_(j) = 1.0 / weights_(j);
            }
        }

        VectorN nodes_;
        VectorN weights_;

        DimensionTpl<N> NDim_;

    }; // class BarycentricInterpolatorTpl

} // namespace galileo

#endif // __galileo_common_math_barycentric_interpolator_hpp__
