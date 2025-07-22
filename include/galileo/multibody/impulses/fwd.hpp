#ifndef __galileo_multibody_impulses_fwd_hpp__
#define __galileo_multibody_impulses_fwd_hpp__

#include "galileo/multibody/fwd.hpp"

namespace galileo
{

    struct ImpulseModelVoid
    {
    }; // struct ImpulseModelVoid

    struct ImpulseDataVoid
    {
    }; // struct ImpulseDataVoid

    template <typename PhaseSpec>
    struct ImpulseModel3dTpl;
    template <typename PhaseSpec>
    struct ImpulseData3dTpl;

    template <typename PhaseSpec>
    struct ImpulseModel6dTpl;
    template <typename PhaseSpec>
    struct ImpulseData6dTpl;

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct ImpulseModelTpl;

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct ImpulseDataTpl;

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    class ImpulseModelManagerTpl;

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    class ImpulseDataManagerTpl;

} // namespace galileo

#endif // __galileo_multibody_impulses_fwd_hpp__
