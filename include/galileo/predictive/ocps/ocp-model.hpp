#ifndef __galileo_predictive_ocps_ocp_model_hpp__
#define __galileo_predictive_ocps_ocp_model_hpp__

#include "galileo/predictive/ocps/ocp-base.hpp"
#include "galileo/predictive/phases/phase-generic.hpp"

namespace galileo
{

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    class OcpModelTpl
    {
    public:
        using BS = BasicSpec;

        using Model_t = OcpModelTpl<BS, PhaseCollectionTpl>;
        using Data_t = OcpDataTpl<BS, PhaseCollectionTpl>;

        using PhaseModel_t = PhaseModelTpl<BS, PhaseCollectionTpl>;
        using PhaseData_t = PhaseDataTpl<BS, PhaseCollectionTpl>;

        OcpModelTpl() {}

        std::vector<PhaseModel_t> phases;

    }; // class OcpModelTpl

} // namespace galileo

#endif // __galileo_predictive_ocps_ocp_model_hpp__
