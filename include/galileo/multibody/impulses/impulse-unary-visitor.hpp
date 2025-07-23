#ifndef __galileo_multibody_impulses_impulse_unary_visitor_hpp__
#define __galileo_multibody_impulses_impulse_unary_visitor_hpp__

#include "galileo/common/meta/unary-visitor.hpp"
#include "galileo/multibody/impulses/impulse-base.hpp"

namespace galileo
{
    namespace fusion
    {
        // Family tag for impulses
        struct ImpulseFamily {};

        // Trait specialization for Impulse family
        template <>
        struct UnaryVisitorFamilyTraits<ImpulseFamily>
        {
            template <typename PS, template <typename> class CollectionTpl>
            using ModelTpl = ImpulseModelTpl<PS, CollectionTpl>;

            template <typename PS, template <typename> class CollectionTpl>
            using DataTpl = ImpulseDataTpl<PS, CollectionTpl>;

            template <typename ModelType, typename PS>
            using ModelBase = ImpulseModelBase<ModelType, PS>;

            template <typename DataType, typename PS>
            using DataBase = ImpulseDataBase<DataType, PS>;
        };

        // Impulse-specific unary visitor base
        template <typename ImpulseVisitorDerived, typename ReturnType = void>
        using ImpulseUnaryVisitorBase = UnaryVisitorBase<ImpulseFamily, ImpulseVisitorDerived, ReturnType>;

    } // namespace fusion

} // namespace galileo

#endif // __galileo_multibody_impulses_impulse_unary_visitor_hpp__
