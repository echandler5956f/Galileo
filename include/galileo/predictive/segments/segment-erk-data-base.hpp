#ifndef __galileo_predictive_segments_segment_erk_data_base_hpp__
#define __galileo_predictive_segments_segment_erk_data_base_hpp__

#include "galileo/predictive/segments/segment-erk-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct SegmentERKDataBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using XNext_t = ArenaMatrixTpl<typename PS::XNext_t>;
        using XNextx_t = ArenaMatrixTpl<typename PS::XNextx_t>;
        using XNextw_t = ArenaMatrixTpl<typename PS::XNextw_t>;
        using L_t = typename PS::L_t;
        using Lx_t = ArenaMatrixTpl<typename PS::Lx_t>;
        using Lw_t = ArenaMatrixTpl<typename PS::Lw_t>;
        using Lxx_t = ArenaMatrixTpl<typename PS::Lxx_t>;
        using Lxw_t = ArenaMatrixTpl<typename PS::Lxw_t>;
        using Lww_t = ArenaMatrixTpl<typename PS::Lww_t>;
        using H_t = ArenaMatrixTpl<typename PS::H_t>;
        using Hx_t = ArenaMatrixTpl<typename PS::Hx_t>;
        using Hw_t = ArenaMatrixTpl<typename PS::Hw_t>;
        using G_t = ArenaMatrixTpl<typename PS::G_t>;
        using Gx_t = ArenaMatrixTpl<typename PS::Gx_t>;
        using Gw_t = ArenaMatrixTpl<typename PS::Gw_t>;

        FORWARD_ACCESSOR(XNext_t, XNext);
        FORWARD_ACCESSOR(XNextx_t, XNextx);
        FORWARD_ACCESSOR(XNextw_t, XNextw);
        FORWARD_ACCESSOR(L_t, L);
        FORWARD_ACCESSOR(Lx_t, Lx);
        FORWARD_ACCESSOR(Lw_t, Lw);
        FORWARD_ACCESSOR(Lxx_t, Lxx);
        FORWARD_ACCESSOR(Lxw_t, Lxw);
        FORWARD_ACCESSOR(Lww_t, Lww);
        FORWARD_ACCESSOR(H_t, H);
        FORWARD_ACCESSOR(Hx_t, Hx);
        FORWARD_ACCESSOR(Hw_t, Hw);
        FORWARD_ACCESSOR(G_t, G);
        FORWARD_ACCESSOR(Gx_t, Gx);
        FORWARD_ACCESSOR(Gw_t, Gw);

    protected:
        inline SegmentERKDataBase() {}
        inline SegmentERKDataBase(const SegmentERKDataBase &clone) { *this = clone; }
        inline SegmentERKDataBase &operator=(const SegmentERKDataBase &clone) { return *this; }

    }; // struct SegmentERKDataBase

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_data_base_hpp__
