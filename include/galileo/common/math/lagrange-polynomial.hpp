#ifndef __galileo_common_math_lagrange_polynomial_hpp__
#define __galileo_common_math_lagrange_polynomial_hpp__

#include "galileo/common/fwd.hpp"

#include <cassert>
#include <iostream>
#include <limits>

namespace galileo
{

    template <typename _NumScalar, int _Options = 0>
    class LagrangePolynomialTpl
    {
    public:
        using NumScalar = _NumScalar;
        static constexpr int Options = _Options;
        using Polynomial = LagrangePolynomialTpl<NumScalar, Options>;
        using MatrixX = Eigen::GMatrix<NumScalar, Eigen::Dynamic, 1, Options>;

        template <typename InputVectorType>
        LagrangePolynomialTpl(const Eigen::MatrixBase<InputVectorType> &coeffs) : coeffs_(coeffs), N_(coeffs.size())
        {
        }

        NumScalar evaluate(const NumScalar &t) const
        {
            // Polynomial evaluation using Horner's method
            // coeffs_[i] corresponds to coefficient of x^(N-1-i)
            // For polynomial: a_0*x^(n-1) + a_1*x^(n-2) + ... + a_(n-1)*x^0
            // Horner's method: ((a_0*x + a_1)*x + a_2)*x + ... + a_(n-1)
            if (N_ == 0)
            {
                return NumScalar(0.0);
            }

            NumScalar result = coeffs_(0);
            for (int i = 1; i < N_; ++i)
            {
                result = result * t + coeffs_(i);
            }

            return result;
        }

        Polynomial derivative() const
        {
            if (N_ <= 1)
            {
                MatrixX coeffs_derivative = MatrixX::Zero(1, 1);
                return Polynomial(coeffs_derivative);
            }

            MatrixX coeffs_derivative = MatrixX::Zero(N_ - 1, 1);
            for (int i = 0; i < N_ - 1; ++i)
            {
                int power = N_ - 1 - i; // power of x for coeffs_[i]
                coeffs_derivative(i) = coeffs_(i) * power;
            }

            return Polynomial(coeffs_derivative);
        }

        Polynomial integral() const
        {
            MatrixX coeffs_integral = MatrixX::Zero(N_ + 1, 1);

            // For descending order coefficients: coeffs_[i] has power (N-1-i)
            // Integration increases power by 1, so new power is (N-i)
            for (int i = 0; i < N_; ++i)
            {
                int power = N_ - 1 - i;    // original power of x for coeffs_[i]
                int new_power = power + 1; // power after integration
                coeffs_integral(i) = coeffs_(i) / new_power;
            }
            coeffs_integral(N_) = 0.0; // constant of integration (lowest power term)

            return Polynomial(coeffs_integral);
        }

        NumScalar integrate(const NumScalar &a, const NumScalar &b) const
        {
            Polynomial indefinite_integral = integral();
            return indefinite_integral.evaluate(b) - indefinite_integral.evaluate(a);
        }

        NumScalar operator()(const NumScalar &t) const { return evaluate(t); }

        template <typename OtherPolynomial>
        Polynomial operator+(const OtherPolynomial &other) const
        {
            int NewN = std::max(N_, other.get_N());
            MatrixX result_coeffs = MatrixX::Zero(NewN, 1);

            MatrixX other_coeffs = other.get_coeffs();

            for (int i = 0; i < N_; ++i)
            {
                result_coeffs(i + NewN - N_) += coeffs_(i);
            }
            for (int i = 0; i < other.get_N(); ++i)
            {
                result_coeffs(i + NewN - other.get_N()) += other_coeffs(i);
            }

            return Polynomial(result_coeffs);
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
            MatrixX result_coeffs = MatrixX::Zero(NewN, 1);

            MatrixX other_coeffs = other.get_coeffs();

            for (int i = 0; i < N_; ++i)
            {
                result_coeffs(i + NewN - N_) += coeffs_(i);
            }
            for (int i = 0; i < other.get_N(); ++i)
            {
                result_coeffs(i + NewN - other.get_N()) -= other_coeffs(i);
            }
            return Polynomial(result_coeffs);
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
            MatrixX result_coeffs = MatrixX::Zero(NewN, 1);

            MatrixX other_coeffs = other.get_coeffs();

            for (int i = 0; i < N_; ++i)
            {
                for (int j = 0; j < other.get_N(); ++j)
                {
                    result_coeffs(i + j) += coeffs_(i) * other_coeffs(j);
                }
            }
            return Polynomial(result_coeffs);
        }

        template <typename OtherPolynomial>
        Polynomial &operator*=(const OtherPolynomial &other)
        {
            *this = operator*(other);
            return *this;
        }

        const MatrixX &get_coeffs() const { return coeffs_; }
        int get_N() const { return N_; }

        friend std::ostream &operator<<(std::ostream &os, const LagrangePolynomialTpl &poly)
        {
            os << "LagrangePolynomial(degree=" << poly.N_ - 1 << ", coefficients=[";
            for (int i = 0; i < poly.N_; ++i)
            {
                if (i > 0) os << ", ";
                os << poly.coeffs_(i);
            }
            os << "])";
            return os;
        }

    protected:
        MatrixX coeffs_;
        int N_;

    }; // class LagrangePolynomialTpl

} // namespace galileo

#endif // __galileo_common_math_lagrange_polynomial_hpp__
