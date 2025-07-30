#ifndef __galileo_predictive_phases_phase_data_base_hpp__
#define __galileo_predictive_phases_phase_data_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"

#define GALILEO_PHASE_DATA_TYPEDEF(Phase)                      \
    using XNext_t = typename traits<Phase>::Data_t::XNext_t;   \
    using XNextx_t = typename traits<Phase>::Data_t::XNextx_t; \
    using XNextw_t = typename traits<Phase>::Data_t::XNextw_t; \
    using L_t = typename traits<Phase>::Data_t::L_t;           \
    using Lx_t = typename traits<Phase>::Data_t::Lx_t;         \
    using Lw_t = typename traits<Phase>::Data_t::Lw_t;         \
    using Lxx_t = typename traits<Phase>::Data_t::Lxx_t;       \
    using Lxw_t = typename traits<Phase>::Data_t::Lxw_t;       \
    using Lww_t = typename traits<Phase>::Data_t::Lww_t;       \
    using H_t = typename traits<Phase>::Data_t::H_t;           \
    using Hx_t = typename traits<Phase>::Data_t::Hx_t;         \
    using Hw_t = typename traits<Phase>::Data_t::Hw_t;         \
    using G_t = typename traits<Phase>::Data_t::G_t;           \
    using Gx_t = typename traits<Phase>::Data_t::Gx_t;         \
    using Gw_t = typename traits<Phase>::Data_t::Gw_t;

namespace galileo
{

    template <typename Derived, typename BasicSpec>
    struct PhaseDataBase
        : public internal::CRTP<Derived>
    {
    public:
        using BS = BasicSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_PHASE_DATA_TYPEDEF(Meta_t);

        XNext_t &XNext_at_seg_i(const int i)
        {
            return this->derived().XNext_at_seg_i_accessor(i);
        }
        const XNext_t &XNext_at_seg_i(const int i) const
        {
            return this->derived().XNext_at_seg_i_accessor(i);
        }

        XNextx_t &XNextx_at_seg_i(const int i)
        {
            return this->derived().XNextx_at_seg_i_accessor(i);
        }
        const XNextx_t &XNextx_at_seg_i(const int i) const
        {
            return this->derived().XNextx_at_seg_i_accessor(i);
        }

        XNextw_t &XNextw_at_seg_i(const int i)
        {
            return this->derived().XNextw_at_seg_i_accessor(i);
        }
        const XNextw_t &XNextw_at_seg_i(const int i) const
        {
            return this->derived().XNextw_at_seg_i_accessor(i);
        }

        L_t &L_at_seg_i(const int i)
        {
            return this->derived().L_at_seg_i_accessor(i);
        }
        const L_t &L_at_seg_i(const int i) const
        {
            return this->derived().L_at_seg_i_accessor(i);
        }

        Lx_t &Lx_at_seg_i(const int i)
        {
            return this->derived().Lx_at_seg_i_accessor(i);
        }
        const Lx_t &Lx_at_seg_i(const int i) const
        {
            return this->derived().Lx_at_seg_i_accessor(i);
        }

        Lw_t &Lw_at_seg_i(const int i)
        {
            return this->derived().Lw_at_seg_i_accessor(i);
        }
        const Lw_t &Lw_at_seg_i(const int i) const
        {
            return this->derived().Lw_at_seg_i_accessor(i);
        }

        Lxx_t &Lxx_at_seg_i(const int i)
        {
            return this->derived().Lxx_at_seg_i_accessor(i);
        }
        const Lxx_t &Lxx_at_seg_i(const int i) const
        {
            return this->derived().Lxx_at_seg_i_accessor(i);
        }

        Lxw_t &Lxw_at_seg_i(const int i)
        {
            return this->derived().Lxw_at_seg_i_accessor(i);
        }
        const Lxw_t &Lxw_at_seg_i(const int i) const
        {
            return this->derived().Lxw_at_seg_i_accessor(i);
        }

        Lww_t &Lww_at_seg_i(const int i)
        {
            return this->derived().Lww_at_seg_i_accessor(i);
        }
        const Lww_t &Lww_at_seg_i(const int i) const
        {
            return this->derived().Lww_at_seg_i_accessor(i);
        }

        H_t &H_at_seg_i(const int i)
        {
            return this->derived().H_at_seg_i_accessor(i);
        }
        const H_t &H_at_seg_i(const int i) const
        {
            return this->derived().H_at_seg_i_accessor(i);
        }

        Hx_t &Hx_at_seg_i(const int i)
        {
            return this->derived().Hx_at_seg_i_accessor(i);
        }
        const Hx_t &Hx_at_seg_i(const int i) const
        {
            return this->derived().Hx_at_seg_i_accessor(i);
        }

        Hw_t &Hw_at_seg_i(const int i)
        {
            return this->derived().Hw_at_seg_i_accessor(i);
        }
        const Hw_t &Hw_at_seg_i(const int i) const
        {
            return this->derived().Hw_at_seg_i_accessor(i);
        }

        G_t &G_at_seg_i(const int i)
        {
            return this->derived().G_at_seg_i_accessor(i);
        }
        const G_t &G_at_seg_i(const int i) const
        {
            return this->derived().G_at_seg_i_accessor(i);
        }

        Gx_t &Gx_at_seg_i(const int i)
        {
            return this->derived().Gx_at_seg_i_accessor(i);
        }
        const Gx_t &Gx_at_seg_i(const int i) const
        {
            return this->derived().Gx_at_seg_i_accessor(i);
        }

        Gw_t &Gw_at_seg_i(const int i)
        {
            return this->derived().Gw_at_seg_i_accessor(i);
        }
        const Gw_t &Gw_at_seg_i(const int i) const
        {
            return this->derived().Gw_at_seg_i_accessor(i);
        }

    protected:
        inline PhaseDataBase(const Model_t &model)
        {
        }

        inline PhaseDataBase(const PhaseDataBase &clone)
        {
        }

        inline PhaseDataBase &operator=(const PhaseDataBase &clone)
        {
            return *this;
        }

    }; // struct PhaseDataBase

} // namespace galileo

#endif // __galileo_predictive_phases_phase_data_base_hpp__
