#ifndef __galileo_predictive_phases_phase_data_base_hpp__
#define __galileo_predictive_phases_phase_data_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct PhaseDataBase
        : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        std::vector<Data_t> &get_segments()
        {
            return segments;
        }
        const std::vector<Data_t> &get_segments() const
        {
            return segments;
        }

        XNext_t &XNext_at_i(const int i)
        {
            return segments[i].XNext;
        }
        const XNext_t &XNext_at_i(const int i) const
        {
            return segments[i].XNext;
        }

        XNextx_t &XNextx_at_i(const int i)
        {
            return segments[i].XNextx;
        }
        const XNextx_t &XNextx_at_i(const int i) const
        {
            return segments[i].XNextx;
        }

        XNextw_t &XNextw_at_i(const int i)
        {
            return segments[i].XNextw;
        }
        const XNextw_t &XNextw_at_i(const int i) const
        {
            return segments[i].XNextw;
        }

        L_t &L_at_i(const int i)
        {
            return segments[i].L;
        }
        const L_t &L_at_i(const int i) const
        {
            return segments[i].L;
        }

        Lx_t &Lx_at_i(const int i)
        {
            return segments[i].Lx;
        }
        const Lx_t &Lx_at_i(const int i) const
        {
            return segments[i].Lx;
        }

        Lw_t &Lw_at_i(const int i)
        {
            return segments[i].Lw;
        }
        const Lw_t &Lw_at_i(const int i) const
        {
            return segments[i].Lw;
        }

        Lxx_t &Lxx_at_i(const int i)
        {
            return segments[i].Lxx;
        }
        const Lxx_t &Lxx_at_i(const int i) const
        {
            return segments[i].Lxx;
        }

        Lxw_t &Lxw_at_i(const int i)
        {
            return segments[i].Lxw;
        }
        const Lxw_t &Lxw_at_i(const int i) const
        {
            return segments[i].Lxw;
        }

        Lww_t &Lww_at_i(const int i)
        {
            return segments[i].Lww;
        }
        const Lww_t &Lww_at_i(const int i) const
        {
            return segments[i].Lww;
        }

        H_t &H_at_i(const int i)
        {
            return segments[i].H;
        }
        const H_t &H_at_i(const int i) const
        {
            return segments[i].H;
        }

        Hx_t &Hx_at_i(const int i)
        {
            return segments[i].Hx;
        }
        const Hx_t &Hx_at_i(const int i) const
        {
            return segments[i].Hx;
        }

        Hw_t &Hw_at_i(const int i)
        {
            return segments[i].Hw;
        }
        const Hw_t &Hw_at_i(const int i) const
        {
            return segments[i].Hw;
        }

        G_t &G_at_i(const int i)
        {
            return segments[i].G;
        }
        const G_t &G_at_i(const int i) const
        {
            return segments[i].G;
        }

        Gx_t &Gx_at_i(const int i)
        {
            return segments[i].Gx;
        }
        const Gx_t &Gx_at_i(const int i) const
        {
            return segments[i].Gx;
        }

        Gw_t &Gw_at_i(const int i)
        {
            return segments[i].Gw;
        }
        const Gw_t &Gw_at_i(const int i) const
        {
            return segments[i].Gw;
        }

        std::vector<Data_t> segments;

    protected:
        inline PhaseDataBase(const Model_t &model)
        {
            segments.reserve(model.get_segments().size());
            for (const auto &segment : model.get_segments())
            {
                segments.push_back(segment.createData());
            }
        }

        inline PhaseDataBase(const PhaseDataBase &clone)
            : segments(clone.segments)
        {
        }

        inline PhaseDataBase &operator=(const PhaseDataBase &clone)
        {
            segments = clone.segments;
            return *this;
        }

    }; // struct PhaseDataBase

} // namespace galileo

#endif // __galileo_predictive_phases_phase_data_base_hpp__
