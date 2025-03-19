#ifndef __galileo_predictive_problem_data_problem_metrics_hpp__
#define __galileo_predictive_problem_data_problem_metrics_hpp__

#include "galileo/predictive/fwd.hpp"

#include "galileo/predictive/problem_data/metrics.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename _VarScalar, typename _NumScalar, int _Options>
        struct ProblemMetrics
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            using Metrics_t = Metrics<VarScalar, NumScalar, Options>;

            std::vector<Metrics_t> running_metrics;
            std::vector<Metrics_t> pre_jump_metrics;
            Metrics_t terminal_metrics;

        }; // class ProblemMetrics

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_problem_data_problem_metrics_hpp__