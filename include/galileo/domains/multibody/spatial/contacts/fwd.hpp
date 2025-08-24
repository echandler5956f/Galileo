#ifndef __galileo_multibody_spatial_contacts_fwd_hpp__
#define __galileo_multibody_spatial_contacts_fwd_hpp__

#include "galileo/domains/multibody/spatial/fwd.hpp"

namespace galileo
{

    struct ContactModelVoid
    {
    }; // struct ContactModelVoid

    struct ContactDataVoid
    {
    }; // struct ContactDataVoid

    template <typename PhaseSpec>
    struct ContactModel3dTpl;
    template <typename PhaseSpec>
    struct ContactData3dTpl;

    template <typename PhaseSpec>
    struct ContactModel6dTpl;
    template <typename PhaseSpec>
    struct ContactData6dTpl;

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct ContactModelTpl;

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct ContactDataTpl;

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    class ContactModelManagerTpl;

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    class ContactDataManagerTpl;

} // namespace galileo

#endif // __galileo_multibody_spatial_contacts_fwd_hpp__
