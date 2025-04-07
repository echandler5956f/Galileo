#ifndef __galileo_core_costs_fwd_hpp__
#define __galileo_core_costs_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    namespace core
    {

        struct CostModelVoid
        {
        }; // struct CostModelVoid

        struct CostDataVoid
        {
        }; // struct CostDataVoid

        template <
            typename PhaseSpec,
            template <typename PS> class ResidualTpl,
            template <typename PS> class ActivationTpl>
        struct CostModelResidualTpl;

        template <
            typename PhaseSpec,
            template <typename PS> class ResidualTpl,
            template <typename PS> class ActivationTpl>
        struct CostDataResidualTpl;

        template <typename PhaseSpec>
        struct CostCollectionDefaultTpl;

        template <
            typename PhaseSpec,
            template <typename PS> class CostCollectionTpl>
        struct CostModelTpl;

        template <
            typename PhaseSpec,
            template <typename PS> class CostCollectionTpl>
        struct CostDataTpl;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_fwd_hpp__