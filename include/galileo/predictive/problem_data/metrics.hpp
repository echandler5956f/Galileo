#ifndef __galileo_predictive_problem_data_metrics_hpp__
#define __galileo_predictive_problem_data_metrics_hpp__

#include "galileo/predictive/fwd.hpp"

namespace galileo
{

    template <typename _VarScalar, typename _NumScalar, int _Option>
    struct Metrics
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using VarScalar = _VarScalar;
        using NumScalar = _NumScalar;
        static constexpr int Options = _Options;

        using VectorXv = Eigen::Matrix<VarScalar, Eigen::Dynamic, Options>;
        using VectorXn = Eigen::Matrix<NumScalar, Eigen::Dynamic, Options>;

        NumScalar cost;
        VectorXn f_feasibility;
        VectorXn h_feasibility;
        VectorXn g_feasibility;

    }; // class Metrics

} // namespace galileo

#endif // __galileo_predictive_problem_data_metrics_hpp__