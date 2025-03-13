#ifndef __galileo_predictive_optimal_control_problem_hpp__
#define __galileo_predictive_optimal_control_problem_hpp__

#include "galileo/predictive/fwd.hpp"

#include "galileo/predictive/trajectory.hpp"

namespace galileo
{

    namespace predictive
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

            using VectorXvs = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1>;
            using VectorXns = Eigen::Matrix<NumScalar, Eigen::Dynamic, 1>;

        protected:
            Trajectory<VarScalar, NumScalar, Options, PhaseCollectionTpl> trajectory_;
            VarScalar cost_;
            VectorXvs x0_;

        }; // class OptimalControlProblem

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_optimal_control_problem_hpp__