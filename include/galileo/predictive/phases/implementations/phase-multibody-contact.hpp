#ifndef __galileo_predictive_phases_phase_multibody_contact_hpp__
#define __galileo_predictive_phases_phase_multibody_contact_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct PhaseMultibodyContactTpl;

    template <typename PhaseSpec>
    struct traits<PhaseMultibodyContactTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseMultibodyContactTpl<PS>;
        using Model_t = PhaseModelMultibodyContactTpl<PS>;
        using Data_t = PhaseDataMultibodyContactTpl<PS>;
    };

    template <typename PhaseSpec>
    struct traits<PhaseDataMultibodyContactTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseMultibodyContactTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct traits<PhaseModelMultibodyContactTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseMultibodyContactTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct PhaseDataMultibodyContactTpl
        : public PhaseDataBase<PhaseDataMultibodyContactTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseMultibodyContactTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseDataBase<PhaseDataMultibodyContactTpl<PS>, PS>;

        PhaseDataMultibodyContactTpl(const Model_t &model)
            : Base(model)
        {
        }

        using Base::get_segments;

        using Base::XNext_at_i;
        using Base::XNextx_at_i;
        using Base::XNextw_at_i;

        using Base::L_at_i;
        using Base::Lx_at_i;
        using Base::Lw_at_i;

        using Base::Lxx_at_i;
        using Base::Lxw_at_i;
        using Base::Lww_at_i;

        using Base::H_at_i;
        using Base::Hx_at_i;
        using Base::Hw_at_i;

        using Base::G_at_i;
        using Base::Gx_at_i;
        using Base::Gw_at_i;

    }; // struct PhaseDataMultibodyContactTpl

    template <typename PhaseSpec>
    class PhaseModelMultibodyContactTpl
        : public PhaseModelBase<PhaseModelMultibodyContactTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseMultibodyContactTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseModelBase<PhaseModelMultibodyContactTpl<PS>, PS>;

        PhaseModelMultibodyContactTpl(const PS &ps)
            : Base(ps)
        {
        }

        using Base::calc;
        using Base::calcDiff;
        using Base::quasiStatic;

        using Base::createData;

        using Base::get_ps;
        using Base::get_segments;

    protected:
    }; // class PhaseModelMultibodyContactTpl

} // namespace galileo

#endif // __galileo_predictive_phases_phase_multibody_contact_hpp__
