#ifndef __galileo_core_data_data_collector_default_hpp__
#define __galileo_core_data_data_collector_default_hpp__

#include "galileo/core/fwd.hpp"
#include <galileo/core/data/data-collector-base.hpp>
#include <galileo/core/data/data-collector-actuation.hpp>
#include <galileo/core/data/data-collector-joint.hpp>
#include <galileo/core/data/data-collector-multibody.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct DataCollectorDefaultTpl : public DataCollectorBase<DataCollectorDefaultTpl<PhaseSpec>>,
                                     public MultibodyDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec>, PhaseSpec>,
                                     public ActuationDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec>, PhaseSpec>,
                                     public JointDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec>, PhaseSpec>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        DataCollectorDefaultTpl(typename PS::RobotData_t *robot,
                                ActuationDataTpl<typename PS::RS> *actuation,
                                JointDataTpl<PS> *joint) : MultibodyDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec>, PhaseSpec>(robot),
                                                           ActuationDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec>, PhaseSpec>(actuation),
                                                           JointDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec>, PhaseSpec>(joint)

        {
        }

    }; // struct DataCollectorDefaultTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_default_hpp__