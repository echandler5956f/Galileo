#ifndef __galileo_common_math_lagrange_polynomial_hpp__
#define __galileo_common_math_lagrange_polynomial_hpp__

#include "galileo/common/fwd.hpp"
#include <cassert>
#include <limits>
#include <iostream>

#include <unsupported/Eigen/Polynomials>

namespace galileo
{
    namespace math
    {

        template <typename _NumScalar, int _Options>
        class LagrangePolynomialTpl
        {
        public:
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;
            using Polynomial = LagrangePolynomialTpl<NumScalar, Options>;
            using GMatrixX = Eigen::GMatrix<NumScalar, Eigen::Dynamic, 1, Options>;

            template <typename InputVectorType>
            LagrangePolynomialTpl(const Eigen::MatrixBase<InputVectorType> &nodes) : nodes_(nodes), N_(nodes.size())
            {
                barycentric_weights_.resize(N_, 1);
                barycentric_weights_.setZero();
            }

            Polynomial derivative() const
            {
                if (N_ <= 1) {
                    GMatrixX nodes_derivative = GMatrixX::Zero(1, 1);
                    return Polynomial(nodes_derivative);
                }

                GMatrixX nodes_derivative = GMatrixX::Zero(N_ - 1, 1);
                for (int i = 0; i < N_ - 1; ++i)
                {
                    int power = N_ - 1 - i;  // power of x for nodes_[i]
                    nodes_derivative(i) = nodes_(i) * power;
                }

                return Polynomial(nodes_derivative);
            }

            Polynomial integral() const
            {
                GMatrixX nodes_integral = GMatrixX::Zero(N_ + 1, 1);

                // For descending order coefficients: nodes_[i] has power (N-1-i)
                // Integration increases power by 1, so new power is (N-i)
                for (int i = 0; i < N_; ++i)
                {
                    int power = N_ - 1 - i;  // original power of x for nodes_[i]
                    int new_power = power + 1;  // power after integration
                    nodes_integral(i) = nodes_(i) / new_power;
                }
                nodes_integral(N_) = 0.0;  // constant of integration (lowest power term)

                return Polynomial(nodes_integral);
            }

            NumScalar integrate(const NumScalar &a, const NumScalar &b) const
            {
                Polynomial indefinite_integral = integral();
                return indefinite_integral.evaluate(b) - indefinite_integral.evaluate(a);
            }

            // Both work, need to decide if we want the Eigen::unsupported dependency
            NumScalar evaluate(const NumScalar &t) const
            {
                return Eigen::poly_eval(nodes_.reverse(), t);

                // Custom polynomial evaluation using Horner's method
                // nodes_[i] corresponds to coefficient of x^(N-1-i)
                // For polynomial: a_0*x^(n-1) + a_1*x^(n-2) + ... + a_(n-1)*x^0
                // Horner's method: ((a_0*x + a_1)*x + a_2)*x + ... + a_(n-1)
                if (N_ == 0) {
                    return NumScalar(0.0);
                }

                NumScalar result = nodes_(0);
                for (int i = 1; i < N_; ++i) {
                    result = result * t + nodes_(i);
                }

                return result;
            }

            NumScalar operator()(const NumScalar &t) const
            {
                return evaluate(t);
            }

            template <typename OtherPolynomial>
            Polynomial operator+(const OtherPolynomial &other) const
            {
                int NewN = std::max(N_, other.get_N());
                GMatrixX result_nodes = GMatrixX::Zero(NewN, 1);

                GMatrixX other_nodes = other.get_nodes();

                for (int i = 0; i < N_; ++i)
                {
                    result_nodes(i + NewN - N_) += nodes_(i);
                }
                for (int i = 0; i < other.get_N(); ++i)
                {
                    result_nodes(i + NewN - other.get_N()) += other_nodes(i);
                }

                return Polynomial(result_nodes);
            }

            template <typename OtherPolynomial>
            Polynomial &operator+=(const OtherPolynomial &other)
            {
                *this = operator+(other);
                return *this;
            }

            template <typename OtherPolynomial>
            Polynomial operator-(const OtherPolynomial &other) const
            {
                int NewN = std::max(N_, other.get_N());
                GMatrixX result_nodes = GMatrixX::Zero(NewN, 1);

                GMatrixX other_nodes = other.get_nodes();

                for (int i = 0; i < N_; ++i)
                {
                    result_nodes(i + NewN - N_) += nodes_(i);
                }
                for (int i = 0; i < other.get_N(); ++i)
                {
                    result_nodes(i + NewN - other.get_N()) -= other_nodes(i);
                }
                return Polynomial(result_nodes);
            }

            template <typename OtherPolynomial>
            Polynomial &operator-=(const OtherPolynomial &other)
            {
                *this = operator-(other);
                return *this;
            }

            template <typename OtherPolynomial>
            Polynomial operator*(const OtherPolynomial &other) const
            {
                int NewN = N_ + other.get_N() - 1;
                GMatrixX result_nodes = GMatrixX::Zero(NewN, 1);

                GMatrixX other_nodes = other.get_nodes();

                for (int i = 0; i < N_; ++i)
                {
                    for (int j = 0; j < other.get_N(); ++j)
                    {
                        result_nodes(i + j) += nodes_(i) * other_nodes(j);
                    }
                }
                return Polynomial(result_nodes);
            }

            template <typename OtherPolynomial>
            Polynomial &operator*=(const OtherPolynomial &other)
            {
                *this = operator*(other);
                return *this;
            }

            // Interpolate a polynomial at t with values w at the nodes, yielding the vector u(t)
            template <typename InputMatrixType, typename OutputVectorType>
            void barycentricInterpolation(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &w, Eigen::MatrixBase<OutputVectorType> &u) const
            {
                assert(w.cols() == N_);
                assert(t >= NumScalar(0.0) && t <= NumScalar(1.0));

                for (int i = 0; i < N_; ++i)
                {
                    if (std::abs(t - nodes_(i)) < std::numeric_limits<NumScalar>::epsilon())
                    {
                        u = w.col(i);
                        return;
                    }
                }

                GMatrixX c = GMatrixX::Zero(N_, 1);
                NumScalar sum_c = 0.0;
                for (int i = 0; i < N_; ++i)
                {
                    c(i) = barycentric_weights_(i) / (t - nodes_(i));
                    sum_c += c(i);
                }

                assert(std::abs(sum_c) > std::numeric_limits<NumScalar>::epsilon() && "Error: Division by zero in BarycentricInterpolation");

                u.setZero();
                for (int i = 0; i < N_; ++i)
                {
                    u += c(i) * w.col(i);
                }
                u /= sum_c;
            }

            template <typename InputMatrixType, typename OutputMatrixType>
            void barycentricInterpolationDiff(const NumScalar &t, const Eigen::MatrixBase<InputMatrixType> &w, Eigen::MatrixBase<OutputMatrixType> &du_dw) const
            {
                assert(w.cols() == N_);
                assert(t >= NumScalar(0.0) && t <= NumScalar(1.0));

                assert(du_dw.rows() == w.rows());
                assert(du_dw.cols() == w.rows() * N_);

                // If t is very close to one of the nodes, the interpolation directly returns w.col(i).
                // In that case, the sensitivity with respect to that column is the identity,
                // and with respect to all other columns is zero.
                for (int i = 0; i < N_; ++i)
                {
                    if (std::abs(t - nodes_[i]) < std::numeric_limits<NumScalar>::epsilon())
                    {
                        du_dw.setZero();
                        du_dw.block(0, i * w.rows(), w.rows(), w.rows()).setIdentity();
                        return;
                    }
                }

                // Compute the barycentric coefficients c_i and their sum.
                GMatrixX c = GMatrixX::Zero(N_, 1);
                NumScalar sum_c = 0.0;
                for (int i = 0; i < N_; ++i)
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
                for (int j = 0; j < N_; ++j)
                {
                    // For each column j, the sensitivity matrix is diagonal with constant c[j] / sum_c.
                    du_dw.block(0, j * w.rows(), w.rows(), w.rows()).diagonal().setConstant(c[j] / sum_c);
                }
            }

            const GMatrixX &get_nodes() const
            {
                return nodes_;
            }

            int get_N() const
            {
                return N_;
            }

            friend std::ostream &operator<<(std::ostream &os, const LagrangePolynomialTpl &poly)
            {
                os << "LagrangePolynomial(degree=" << poly.N_ - 1 << ", coefficients=[";
                for (int i = 0; i < poly.N_; ++i)
                {
                    if (i > 0) os << ", ";
                    os << poly.nodes_(i);
                }
                os << "])";
                return os;
            }

        protected:
            void compute_barycentric_weights()
            {
                /*Barycentric weights*/
                for (int j = 0; j < N_; ++j)
                    barycentric_weights_(j) = 1.0;

                /*For all collocation points*/
                for (int j = 0; j < N_; ++j)
                {
                    for (int r = 0; r < N_; ++r)
                    {
                        if (r != j)
                        {
                            barycentric_weights_(j) *= (nodes_(j) - nodes_(r));
                        }
                    }
                    barycentric_weights_(j) = 1.0 / barycentric_weights_(j);
                }
            }

            GMatrixX nodes_;
            GMatrixX barycentric_weights_;
            int N_;

        }; // class LagrangePolynomialTpl

    } // namespace math

} // namespace galileo

#endif // __galileo_common_math_lagrange_polynomial_hpp__
