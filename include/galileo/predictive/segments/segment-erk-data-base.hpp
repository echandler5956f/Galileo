#ifndef __galileo_predictive_segments_segment_erk_data_base_hpp__
#define __galileo_predictive_segments_segment_erk_data_base_hpp__

#include "galileo/predictive/segments/segment-erk-base.hpp"
#include "galileo/predictive/segments/segment-erk-model-base.hpp"

namespace galileo
{
    namespace predictive
    {

        template <typename Derived, typename PhaseSpec>
        struct SegmentERKDataBase : internal::CRTP<SegmentERKDataBase<Derived, PhaseSpec>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

            using SegmentERKDerived = typename traits<Derived>::SegmentERKDerived;
            using SegmentERKDataDerived = typename traits<SegmentERKDerived>::SegmentERKDataDerived;
            using SegmentERKModelDerived = typename traits<SegmentERKDerived>::SegmentERKModelDerived;

            FORWARD_ACCESSOR(XNext_t, XNext);
            FORWARD_ACCESSOR(Fx_t, Fx);
            FORWARD_ACCESSOR(Fw_t, Fw);

            FORWARD_ACCESSOR(L_t, L);
            FORWARD_ACCESSOR(Lx_t, Lx);
            FORWARD_ACCESSOR(Lw_t, Lw);
            FORWARD_ACCESSOR(Lxx_t, Lxx);
            FORWARD_ACCESSOR(Lxw_t, Lxw);
            FORWARD_ACCESSOR(Lww_t, Lww);

            FORWARD_ACCESSOR(H_t, H);
            FORWARD_ACCESSOR(Hx_t, Hx);
            FORWARD_ACCESSOR(Hu_t, Hu);
            FORWARD_ACCESSOR(G_t, G);
            FORWARD_ACCESSOR(Gx_t, Gx);
            FORWARD_ACCESSOR(Gu_t, Gu);

        protected:
            inline SegmentERKDataBase()
            {
            }

            inline SegmentERKDataBase(const SegmentERKDataBase &clone)
            {
                *this = clone;
            }

            inline SegmentERKDataBase &operator=(const SegmentERKDataBase &clone)
            {
                return *this;
            }

        }; // struct SegmentERKDataBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_data_base_hpp__
