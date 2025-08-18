#ifndef __galileo_predictive_phases_phase_visitor_base_hpp__
#define __galileo_predictive_phases_phase_visitor_base_hpp__

#include "galileo/common/visitors/unary-visitor.hpp"
#include "galileo/common/visitors/binary-visitor.hpp"
#include "galileo/predictive/phases/phase-base.hpp"

namespace galileo
{
    namespace fusion
    {
        // Family tag for phases
        struct PhaseFamily
        {
        };

        // Trait specialization for Phase family
        template <>
        struct UnaryVisitorFamilyTraits<PhaseFamily>
        {
            template <typename BS, template <typename> class CollectionTpl>
            using ModelTpl = PhaseModelTpl<BS, CollectionTpl>;

            template <typename BS, template <typename> class CollectionTpl>
            using DataTpl = PhaseDataTpl<BS, CollectionTpl>;

            template <typename ModelType, typename BS>
            using ModelBase = PhaseModelBase<ModelType, BS>;

            template <typename DataType, typename BS>
            using DataBase = PhaseDataBase<DataType, BS>;
        };

        // Phase-specific unary visitor base
        template <typename PhaseVisitorDerived, typename ReturnType = void>
        using PhaseUnaryVisitorBase = UnaryVisitorBase<PhaseFamily, PhaseVisitorDerived, ReturnType>;

        template <typename PhaseVisitorDerived, typename ReturnType = void>
        using PhaseBinaryVisitorBase = BinaryVisitorBase<PhaseFamily, PhaseVisitorDerived, ReturnType>;

    } // namespace fusion

} // namespace galileo

#endif // __galileo_predictive_phases_phase_visitor_base_hpp__
