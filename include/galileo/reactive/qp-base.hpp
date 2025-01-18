#ifndef __galileo_reactive_qp_base_hpp__
#define __galileo_reactive_qp_base_hpp__

#include "galileo/reactive/fwd.hpp"

namespace galileo
{
    namespace reactive
    {
    template <typename _Scalar>
    class QPBase
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using Scalar = _Scalar;
        using VectorXs = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        using MatrixXs = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

        QPBase() = default;

        virtual void solve(const MatrixXs &H, const VectorXs &g, const MatrixXs &A, const VectorXs &lb, const VectorXs &ub, VectorXs &w) = 0;
    };

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_qp_base_hpp__