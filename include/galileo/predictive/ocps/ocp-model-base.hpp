#ifndef __galileo_predictive_ocps_ocp_model_base_hpp__
#define __galileo_predictive_ocps_ocp_model_base_hpp__

#include "galileo/predictive/ocps/ocp-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseCollectionTpl>
    class OcpModelBase
        : public internal::CRTP<Derived>
    {
    public:
        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

    protected:
        inline OcpModelBase()
        {
        }

        inline OcpModelBase(const OcpModelBase &clone)
        {
            *this = clone;
        }

        inline OcpModelBase &operator=(const OcpModelBase &clone)
        {
            return *this;
        }

    }; // class OcpModelBase

} // namespace galileo

#endif // __galileo_predictive_ocps_ocp_model_base_hpp__
