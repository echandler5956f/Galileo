#ifndef __galileo_predictive_nodes_fwd_hpp__
#define __galileo_predictive_nodes_fwd_hpp__

#include "galileo/predictive/fwd.hpp"

namespace galileo
{

    namespace predictive
    {

        template <
            typename VarScalar,
            typename NumScalar,
            int Options,
            int NX,
            int NU,
            int NDX,
            template <typename V, typename N, int O, int NX, int NU, int NDX> class StateTpl,
            template <typename V, typename N, int O> class ActuationTpl,
            template <typename V, typename N, int O> class ConstraintCollectionTpl,
            template <typename V, typename N, int O> class CostCollectionTpl>
        class NodeBaseTpl;

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_nodes_fwd_hpp__