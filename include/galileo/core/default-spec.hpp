#ifndef __galileo_core_default_spec_hpp__
#define __galileo_core_default_spec_hpp__

#include "galileo/core/system-spec.hpp"

namespace galileo
{
    /* ---------------------------------------------------------------- */
    /* Defines the system dimensions and types used in the core library. */
    /* ---------------------------------------------------------------- */
    template <typename BasicSpec,
              int _NQ,
              int _NV,
              int _NUa,
              template <typename> class StateTpl,
              template <typename> class ActuationTpl>
    struct DefaultSpecTpl : public SystemSpecTpl<BasicSpec, _NQ, _NV, _NUa>
    {
        using DS = DefaultSpecTpl<BasicSpec, _NQ, _NV, _NUa, StateTpl, ActuationTpl>;
        using Base = SystemSpecTpl<BasicSpec, _NQ, _NV, _NUa>;

        using BS = BasicSpec;
        GALILEO_SYSTEM_SPEC_MASTER_TYPEDEF(Base);

        /* ---------------------------------------------------------------- */
        /* Template types */
        /* ---------------------------------------------------------------- */
        using State_t = StateTpl<DS>;

        using ActuationMeta_t = ActuationTpl<DS>;
        using ActuationModel_t = typename traits<ActuationMeta_t>::Model_t;
        using ActuationData_t = typename traits<ActuationMeta_t>::Data_t;

        using Base::nq_dim_;
        using Base::nv_dim_;
        using Base::nx_dim_;
        using Base::ndx_dim_;
        using Base::nua_dim_;

        DefaultSpecTpl() : Base() {}

        using Base::get_nq;
        using Base::get_nq_dim;
        using Base::get_nv;
        using Base::get_nv_dim;
        using Base::get_nx;
        using Base::get_nx_dim;
        using Base::get_ndx;
        using Base::get_ndx_dim;
        using Base::get_nua;
        using Base::get_nua_dim;
        using Base::is_valid_spec;
        using Base::display;

        friend std::ostream &operator<<(std::ostream &os, const DefaultSpecTpl &ds)
        {
            os << "DefaultSpec: {\n";
            ds.display(os, "  ");
            os << "}";
            return os;
        }
    };

} // namespace galileo

#endif // __galileo_core_default_spec_hpp__
