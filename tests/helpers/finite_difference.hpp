#ifndef __galileo_tests_helpers_finite_difference_hpp__
#define __galileo_tests_helpers_finite_difference_hpp__

#include <functional>
#include <Eigen/Dense>

namespace galileo
{
    namespace tests
    {

        template <typename Func, typename TangentVec, typename DomainVec, typename RangeVec>
        void compute_finite_difference_jacobian(
            const Func &func,
            const DomainVec &x,
            double epsilon,
            Eigen::MatrixBase<Eigen::Matrix<typename TangentVec::Scalar, RangeVec::RowsAtCompileTime, TangentVec::RowsAtCompileTime>> &J)
        {

            using Scalar = typename TangentVec::Scalar;
            TangentVec perturbation = TangentVec::Zero(x.size());
            RangeVec y_plus, y_minus;

            for (int i = 0; i < x.size(); ++i)
            {
                perturbation(i) = epsilon;
                func(x + perturbation, y_plus);
                func(x - perturbation, y_minus);
                J.col(i) = (y_plus - y_minus) / (2 * epsilon);
                perturbation(i) = 0;
            }
        }

        template <typename Manifold, typename Func, typename DomainVec, typename RangeVec, typename JacobianMat>
        void compute_finite_difference_jacobian_manifold(
            const Manifold &manifold,
            const Func &func, // func(x, y_out)
            const DomainVec &x,
            double epsilon,
            Eigen::MatrixBase<JacobianMat> &J)
        {

            using TangentVec = Eigen::Matrix<typename DomainVec::Scalar, Eigen::Dynamic, 1>;
            TangentVec perturbation = TangentVec::Zero(manifold.get_ndx());
            RangeVec y_plus, y_minus;
            DomainVec x_plus(manifold.get_nx()), x_minus(manifold.get_nx());

            for (int i = 0; i < manifold.get_ndx(); ++i)
            {
                perturbation(i) = epsilon;

                manifold.integrate(x, perturbation, x_plus);
                func(x_plus, y_plus);

                manifold.integrate(x, -perturbation, x_minus);
                func(x_minus, y_minus);

                J.col(i) = (y_plus - y_minus) / (2 * epsilon);
                perturbation(i) = 0;
            }
        }

    } // namespace tests
    
} // namespace galileo

#endif // __galileo_tests_helpers_finite_difference_hpp__