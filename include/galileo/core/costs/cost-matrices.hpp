#ifndef __galileo_core_costs_cost_matrices_hpp__
#define __galileo_core_costs_cost_matrices_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename Scalar, bool ShareData, int NX, int NU, int Options>
        struct CostMatrices
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            // Underlying Eigen matrix types
            using Lx_Mat_t = Eigen::Matrix<Scalar, 1, NX, Options>;
            using Lu_Mat_t = Eigen::Matrix<Scalar, 1, NU, Options>;
            using Lxx_Mat_t = Eigen::Matrix<Scalar, NX, NX, Options>;
            using Lxu_Mat_t = Eigen::Matrix<Scalar, NX, NU, Options>;
            using Luu_Mat_t = Eigen::Matrix<Scalar, NU, NU, Options>;

            using L_t = Scalar;
            using Lx_t = SelectMatrix<Lx_Mat_t, ShareData>;
            using Lu_t = SelectMatrix<Lu_Mat_t, ShareData>;
            using Lxx_t = SelectMatrix<Lxx_Mat_t, ShareData>;
            using Lxu_t = SelectMatrix<Lxu_Mat_t, ShareData>;
            using Luu_t = SelectMatrix<Luu_Mat_t, ShareData>;

            L_t L;
            Lx_t Lx;
            Lu_t Lu;
            Lxx_t Lxx;
            Lxu_t Lxu;
            Luu_t Luu;

        }; // struct CostMatrices

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_matrices_hpp__