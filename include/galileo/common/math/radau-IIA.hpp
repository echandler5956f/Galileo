#ifndef __galileo_common_math_radau_IIA_hpp__
#define __galileo_common_math_radau_IIA_hpp__

#include "galileo/common/fwd.hpp"
#include "galileo/common/math/jacobi-roots.hpp"
#include "galileo/common/math/lagrange-polynomial.hpp"

#include <cmath>
#include <vector>

namespace galileo
{

    constexpr inline size_t binom(size_t n, size_t k) noexcept
    {
        return (k > n) ? 0 : // out of range
                   (k == 0 || k == n) ? 1
                                      : // edge
                   (k == 1 || k == n - 1) ? n
                                          :     // first
                   binom(n - 1, k - 1) * n / k; // recursive
    }

    template <typename _NumScalar, int _N, int _Options = 0>
    class RadauIIATpl
    {
    public:
        using NumScalar = _NumScalar;
        static constexpr int N = _N;
        static constexpr int Options = _Options;

        using Polynomial = LagrangePolynomialTpl<NumScalar, Options>;

        RadauIIATpl() : jacobi_roots_(1.0, 0.0)
        {
        }

        RadauIIATpl(const JacobiRootsTpl<NumScalar, N - 1, Options> &jacobi_roots) : jacobi_roots_(jacobi_roots)
        {
            GALILEO_ASSERT(jacobi_roots_.get_alpha() == 1.0 && jacobi_roots_.get_beta() == 0.0,
                           "RadauIIATpl: Jacobi roots must be from shifted Legendre polynomials (alpha = 1.0, beta = 0.0).");
        }

        void compute_terms()
        {
            jacobi_roots_.compute_roots();
            Eigen::GMatrix<NumScalar, N - 1, 1, Options> roots = jacobi_roots_.get_roots();
            nodes_ << roots, 1.0;

            butcher_matrix_.setZero();
            weights_.setZero();
            W_.setZero();
            X_.setZero();

            polynomials_.clear();

            Eigen::GMatrix<NumScalar, 2, 1, Options> tmp;
            for (int j = 0; j < N; ++j)
            {
                Eigen::GMatrix<NumScalar, 1, 1, Options> scalar_one(1.0);
                Polynomial p = Polynomial(scalar_one);
                for (int r = 0; r < N; ++r)
                {
                    if (r != j)
                    {
                        NumScalar denom = nodes_(j) - nodes_(r);
                        tmp(0) = 1.0 / denom;
                        tmp(1) = -nodes_(r) / denom;
                        p = p * Polynomial(tmp);
                    }
                }
                polynomials_.push_back(p);
            }

            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    Polynomial p_j = polynomials_[j];
                    butcher_matrix_(i, j) = p_j.integrate(0.0, nodes_(i));
                    W_(i, j) = shifted_legendre_ortho(j, nodes_(i));
                }
            }

            weights_ = butcher_matrix_.template bottomRows<1>().transpose();

            if (N > 0)
            {
                X_(0, 0) = 0.5;
                if (N > 1)
                {
                    for (int i = 0; i < N - 1; ++i)
                    {
                        NumScalar zeta_k = calc_zeta(i);
                        X_(i, i + 1) = -zeta_k;
                        X_(i + 1, i) = zeta_k;
                    }
                    X_.template bottomRightCorner<1, 1>()(0) = 1.0 / ((4.0 * N) - 2.0);
                }
            }
        }

        NumScalar shifted_legendre_ortho(int k, NumScalar x) const
        {
            GALILEO_ASSERT(k >= 0, "RadauIIATpl: Degree k for shifted Legendre polynomial must be non-negative.");

            if (k == 0)
            {
                return 1.0;
            }
            NumScalar term_sum = 0.0;
            for (int j = 0; j <= k; ++j)
            {
                NumScalar comb_k_j = binom(k, j);
                NumScalar comb_jpk_j = binom(j + k, j);
                NumScalar term = std::pow(-1.0, j + k) * comb_k_j * comb_jpk_j * std::pow(x, j);
                term_sum += term;
            }
            NumScalar normalization_factor = std::sqrt(2.0 * k + 1.0);
            return normalization_factor * term_sum;
        }

        NumScalar calc_zeta(int i) const
        {
            return 1.0 / (2.0 * std::sqrt(((2.0 * i) + 1.0) * ((2.0 * i) + 3.0)));
        }

        const Eigen::GMatrix<NumScalar, N, 1, Options> &get_nodes() const
        {
            return nodes_;
        }

        const Eigen::GMatrix<NumScalar, N, N, Options> &get_butcher_matrix() const
        {
            return butcher_matrix_;
        }

        const Eigen::GMatrix<NumScalar, N, 1, Options> &get_weights() const
        {
            return weights_;
        }

        const Eigen::GMatrix<NumScalar, N, N, Options> &get_W() const
        {
            return W_;
        }

        const Eigen::GMatrix<NumScalar, N, N, Options> &get_X() const
        {
            return X_;
        }

        const JacobiRootsTpl<NumScalar, N - 1, Options> &get_jacobi_roots() const
        {
            return jacobi_roots_;
        }

    protected:
        JacobiRootsTpl<NumScalar, N - 1, Options> jacobi_roots_;
        Eigen::GMatrix<NumScalar, N, N, Options> butcher_matrix_;
        Eigen::GMatrix<NumScalar, N, 1, Options> weights_;
        Eigen::GMatrix<NumScalar, N, 1, Options> nodes_;
        Eigen::GMatrix<NumScalar, N, N, Options> W_;
        Eigen::GMatrix<NumScalar, N, N, Options> X_;
        std::vector<Polynomial> polynomials_;

    }; // class LagrangePolynomialTpl

} // namespace galileo

#endif // __galileo_common_math_radau_IIA_hpp__
