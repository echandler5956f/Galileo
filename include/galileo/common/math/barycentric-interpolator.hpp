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
        BarycentricInterpolatorTpl(const Eigen::MatrixBase<InputVectorType> &nodes) : nodes_(nodes)
        {
            NDim_.set_value(nodes_.size());
            compute_weights();
        }

        template <typename InputMatrixType, typename OutputVectorType>
        void calc(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &w, Eigen::MatrixBase<OutputVectorType> &u) const
        {
            using RowVectorType = typename Eigen::internal::plain_row_type<OutputVectorType>::type;

            GALILEO_ASSERT(w.cols() == NDim_.value(), "BarycentricInterpolator: Input matrix must have the same number of columns as the number of nodes.");
            GALILEO_ASSERT(t >= NumScalar(0.) && t <= NumScalar(1.), "BarycentricInterpolator: Time must be between 0 and 1.");

            // Compute the interpolated value
            RowVectorType numerator = RowVectorType::Zero(NDim_.value());
            RowVectorType denominator = RowVectorType::Zero(NDim_.value());
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
                denominator += RowVectorType::Constant(NDim_.value(), interpolant);
            }

            GALILEO_ASSERT((denominator.array() == 0).any() == false, "BarycentricInterpolator: Denominator is zero.");

            u = numerator.array() / denominator.array();
        }

        template <typename InputMatrixType, typename OutputMatrixType>
        void calcDiff(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &w, Eigen::MatrixBase<OutputMatrixType> &du_dw) const
        {
            GALILEO_ASSERT(w.cols() == NDim_.value(), "BarycentricInterpolator: Input matrix must have the same number of columns as the number of nodes.");
            GALILEO_ASSERT(t >= NumScalar(0.) && t <= NumScalar(1.), "BarycentricInterpolator: Time must be between 0 and 1.");

            GALILEO_ASSERT(du_dw.rows() == NDim_.value(), "BarycentricInterpolator: Output matrix must have the same number of rows as the number of nodes.");
            GALILEO_ASSERT(du_dw.cols() == NDim_.value() * NDim_.value(), "BarycentricInterpolator: Output matrix must have the same number of columns as the number of nodes squared.");

            // If t is very close to one of the nodes, the interpolation directly returns w.col(i).
            // In that case, the sensitivity with respect to that column is the identity,
            // and with respect to all other columns is zero.
            for (int i = 0; i < NDim_.value(); ++i)
            {
                if (std::abs(t - nodes_(i)) < std::numeric_limits<NumScalar>::epsilon())
                {
                    du_dw.setZero();
                    block(du_dw, 0, i * NDim_.value(), NDim_, NDim_) = MatrixN::Identity(NDim_.value(), NDim_.value());
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
            GALILEO_ASSERT(std::abs(sum_c) > std::numeric_limits<NumScalar>::epsilon(), "BarycentricInterpolator: Sum of coefficients is zero.");

            // The interpolated value is u = (sum_i c_i * w.col(i)) / sum_c.
            // Thus, for each element k of u:
            //     u(k) = (sum_i c_i * w(k,i)) / sum_c.
            // Therefore, the partial derivative with respect to w(k,j) is c_j/sum_c (if the row index matches).
            // We pack these derivatives into a vector of matrices, where each matrix is (w.rows() x w.rows())
            // representing the derivative with respect to one column of w.
            for (int j = 0; j < NDim_.value(); ++j)
            {
                // For each column j, the sensitivity matrix is diagonal with constant c[j] / sum_c.
                block(du_dw, 0, j * NDim_.value(), NDim_, NDim_).diagonal().setConstant(c(j) / sum_c);
            }
        }

        constexpr int get_n() const
        {
            return NDim_.value();
        }

        const DimensionTpl<N> &NDim() const
        {
            return NDim_;
        }

    protected:
        void compute_weights()
        {
            weights_ = VectorN::Zero(NDim_.value());

            /*Barycentric weights*/
            for (int j = 0; j < NDim_.value(); ++j)
                weights_(j) = 1.0;

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
