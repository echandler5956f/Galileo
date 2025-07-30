#ifndef __galileo_predictive_ocps_ocp_data_base_hpp__
#define __galileo_predictive_ocps_ocp_data_base_hpp__

#include "galileo/predictive/ocps/ocp-base.hpp"

namespace galileo
{

    template <typename Derived>
    struct OcpDataBase
        : public internal::CRTP<Derived>
    {
    public:
        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

    protected:
        inline OcpDataBase(const Model_t &model)
        {
        }

        inline OcpDataBase(const OcpDataBase &clone)
        {
        }

        inline OcpDataBase &operator=(const OcpDataBase &clone)
        {
            return *this;
        }

    }; // struct OcpDataBase

} // namespace galileo

#endif // __galileo_predictive_ocps_ocp_data_base_hpp__
