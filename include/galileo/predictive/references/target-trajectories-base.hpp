#ifndef __galileo_predictive_references_target_trajectories_base_hpp__
#define __galileo_predictive_references_target_trajectories_base_hpp__

#include "galileo/predictive/fwd.hpp"

#define GALILEO_TARGET_TRAJECTORIES_BASIC_TYPEDEF(TargetTrajectories) \
    using VarScalar = typename traits<TargetTrajectories>::VarScalar;    \
    using NumScalar = typename traits<TargetTrajectories>::NumScalar;    \
    static constexpr int Options = traits<TargetTrajectories>::Options;

#define GALILEO_TARGET_TRAJECTORIES_TYPEDEF(TargetTrajectories)    \
    using VectorXv = typename traits<TargetTrajectories>::VectorXv; \
    using VectorXn = typename traits<TargetTrajectories>::VectorXn;

namespace galileo
{

    namespace predictive
    {

        template <typename Derived>
        class TargetTrajectoriesBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using TargetTrajectoriesDerived = typename traits<Derived>::TargetTrajectoriesDerived;
            GALILEO_TARGET_TRAJECTORIES_BASIC_TYPEDEF(TargetTrajectoriesDerived);
            GALILEO_TARGET_TRAJECTORIES_TYPEDEF(TargetTrajectoriesDerived);

            VectorXv getDesiredState(NumScalar t) const
            {
                return derived().getDesiredState(t);
            }

            VectorXn getDesiredInput(NumScalar t) const
            {
                return derived().getDesiredInput(t);
            }

        }; // class TargetTrajectoriesBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_references_target_trajectories_base_hpp__
