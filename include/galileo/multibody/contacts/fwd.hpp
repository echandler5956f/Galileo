#ifndef __galileo_multibody_contacts_fwd_hpp__
#define __galileo_multibody_contacts_fwd_hpp__

#include "galileo/multibody/fwd.hpp"

namespace galileo
{

    struct ContactModelVoid
    {
    }; // struct ContactModelVoid

    struct ContactDataVoid
    {
    }; // struct ContactDataVoid

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactModelTpl;

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactDataTpl;

} // namespace galileo

#endif // __galileo_multibody_contacts_fwd_hpp__