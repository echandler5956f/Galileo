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
        struct PhaseFamily {};

        // Trait specialization for Phase family
        template <>
        struct UnaryVisitorFamilyTraits<PhaseFamily>
        {
            template <typename PS, template <typename> class CollectionTpl>
            using ModelTpl = PhaseModelTpl<PS, CollectionTpl>;

            template <typename PS, template <typename> class CollectionTpl>
            using DataTpl = PhaseDataTpl<PS, CollectionTpl>;

            template <typename ModelType, typename PS>
            using ModelBase = PhaseModelBase<ModelType, PS>;

            template <typename DataType, typename PS>
            using DataBase = PhaseDataBase<DataType, PS>;
        };

        // Phase-specific unary visitor base
        template <typename PhaseVisitorDerived, typename ReturnType = void>
        using PhaseUnaryVisitorBase = UnaryVisitorBase<PhaseFamily, PhaseVisitorDerived, ReturnType>;

        template <typename PhaseVisitorDerived, typename ReturnType = void>
        using PhaseBinaryVisitorBase = BinaryVisitorBase<PhaseFamily, PhaseVisitorDerived, ReturnType>;

    } // namespace fusion

} // namespace galileo

#endif // __galileo_predictive_phases_phase_visitor_base_hpp__
