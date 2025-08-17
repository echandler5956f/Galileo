#ifndef __galileo_predictive_ocps_ocp_data_hpp__
#define __galileo_predictive_ocps_ocp_data_hpp__

#include "galileo/predictive/ocps/ocp-base.hpp"

namespace galileo
{

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct OcpDataTpl
    {
    public:
        using BS = BasicSpec;

        using Model_t = OcpModelTpl<BS, PhaseCollectionTpl>;
        using Data_t = OcpDataTpl<BS, PhaseCollectionTpl>;

        using PhaseModel_t = PhaseModelTpl<BS, PhaseCollectionTpl>;
        using PhaseData_t = PhaseDataTpl<BS, PhaseCollectionTpl>;

        OcpDataTpl() {}

        std::vector<PhaseData_t> phases;

    }; // struct OcpDataTpl

} // namespace galileo

#endif // __galileo_predictive_ocps_ocp_data_hpp__
