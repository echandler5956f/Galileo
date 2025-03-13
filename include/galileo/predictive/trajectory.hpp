#ifndef __galileo_predictive_trajectory_hpp__
#define __galileo_predictive_trajectory_hpp__

#include "galileo/predictive/fwd.hpp"
// #include "galileo/predictive/phases/phase-generic.hpp"
#include "galileo/utils/aligned-vector.hpp"

namespace galileo
{

    namespace predictive
    {
        
        template <typename _VarScalar, typename _NumScalar, int _Options, template <typename, typename, int> class PhaseCollectionTpl>
        class Trajectory
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;
            using PhaseCollection = PhaseCollectionTpl<VarScalar, NumScalar, Options>;

            using PhaseModel = PhaseModelTpl<VarScalar, NumScalar, Options>;
            using PhaseData = PhaseDataTpl<VarScalar, NumScalar, Options>;

            using PhaseModelVector = typename GALILEO_ALIGNED_STD_VECTOR(PhaseModel);
            using PhaseDataVector = typename GALILEO_ALIGNED_STD_VECTOR(PhaseData);
        
            protected:
                PhaseModelVector phase_models_;
                PhaseDataVector phase_datas_;

                // Jump map handled via visitor pattern


        }; // class Trajectory

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_trajectory_hpp__