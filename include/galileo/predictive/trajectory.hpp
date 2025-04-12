#ifndef __galileo_predictive_trajectory_hpp__
#define __galileo_predictive_trajectory_hpp__

#include "galileo/predictive/fwd.hpp"
// #include "galileo/predictive/phases/phase-generic.hpp"
#include "galileo/utils/aligned-vector.hpp"

namespace galileo
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

        using PhaseModel = PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>;
        using PhaseData = PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>;

        using PhaseModelVector = typename GALILEO_ALIGNED_STD_VECTOR(PhaseModel);
        using PhaseDataVector = typename GALILEO_ALIGNED_STD_VECTOR(PhaseData);

    protected:
        void buildOffsets()
        {
            std::vector<std::size_t> offset(phase_models_.size() + 1, 0);
            // offset[0] = 0 by default

            for (std::size_t i = 0; i < phase_models_.size(); i++)
            {
                // how many segments in this phase?
                std::size_t count = 0;
                boost::apply_visitor([&count](auto const &p)
                                     { count = p.segments_.size(); }, phase_models_[i]);

                offset[i + 1] = offset[i] + count;
            }
            // offset.back() = total number of segments
            phase_offsets_ = offset;
        }

        size_t find_phase_index(size_t g)
        {
            // offsets[i] <= g < offsets[i+1], find i
            // e.g. by lower_bound:
            //    we want the largest i s.t. offsets[i] <= g
            //    but we must also check g < offsets[i+1].
            // Example:
            auto it = std::upper_bound(phase_offsets_.begin(), phase_offsets_.end(), g);
            // upper_bound returns first element > g, so the phase is (it-1).
            // But we must be sure it's not at offsets.begin()
            size_t i = (it - phase_offsets_.begin()) - 1;
            return i;
        }

        PhaseModelVector phase_models_;
        PhaseDataVector phase_datas_;
        std::vector<std::size_t> phase_offsets_;

        // Jump map handled via visitor pattern

    }; // class Trajectory

} // namespace galileo

#endif // __galileo_predictive_trajectory_hpp__