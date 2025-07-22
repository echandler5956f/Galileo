#ifndef __galileo_multibody_impulses_impulse_collection_hpp__
#define __galileo_multibody_impulses_impulse_collection_hpp__

#include "galileo/multibody/impulses/fwd.hpp"
#include "galileo/multibody/impulses/implementations/impulse-3d.hpp"
#include "galileo/multibody/impulses/implementations/impulse-6d.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct ImpulseCollectionDefaultTpl
    {
        using PS = PhaseSpec;

        using ImpulseModelVariant_t = boost::variant<ImpulseModel3dTpl<PS>, ImpulseModel6dTpl<PS>>;
        using ImpulseDataVariant_t = boost::variant<ImpulseData3dTpl<PS>, ImpulseData6dTpl<PS>>;

    }; // struct ImpulseCollectionDefaultTpl

    template <typename PhaseSpec>
    using ImpulseModelVariantTpl = ImpulseCollectionDefaultTpl<PhaseSpec>::ImpulseModelVariant_t;

    template <typename PhaseSpec>
    using ImpulseDataVariantTpl = ImpulseCollectionDefaultTpl<PhaseSpec>::ImpulseDataVariant_t;

} // namespace galileo

#endif // __galileo_multibody_impulses_impulse_collection_hpp__
