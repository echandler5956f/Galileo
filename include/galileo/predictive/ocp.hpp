#ifndef __galileo_predictive_ocp_hpp__
#define __galileo_predictive_ocp_hpp__

#include "galileo/predictive/fwd.hpp"

#include "galileo/predictive/trajectory.hpp"

namespace galileo
{

    template <typename _VarScalar, typename _NumScalar, int _Options, template <typename, typename, int> class PhaseCollectionTpl>
    class OptimalControlProblem
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using VarScalar = _VarScalar;
        using NumScalar = _NumScalar;
        static constexpr int Options = _Options;
        using PhaseCollection = PhaseCollectionTpl<VarScalar, NumScalar, Options>;

        using VectorXv = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1>;
        using VectorXn = Eigen::Matrix<NumScalar, Eigen::Dynamic, 1>;

        // protected:
        Trajectory<VarScalar, NumScalar, Options, PhaseCollectionTpl> trajectory_;
        VarScalar cost_;
        VectorXv x0_;

    }; // class OptimalControlProblem

} // namespace galileo

#endif // __galileo_predictive_ocp_hpp__