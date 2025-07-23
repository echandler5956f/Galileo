#ifndef __galileo_core_data_data_collector_base_hpp__
#define __galileo_core_data_data_collector_base_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    template <typename Derived>
    struct DataCollectorBase
        : public internal::CRTP<Derived>
    {
    protected:
        inline DataCollectorBase()
        {
        }

        inline DataCollectorBase(const DataCollectorBase &clone)
        {
            *this = clone;
        }

        inline DataCollectorBase &operator=(const DataCollectorBase &clone)
        {
            return *this;
        }

    }; // struct DataCollectorBase

} // namespace galileo

#endif // __galileo_core_data_data_collector_base_hpp__
