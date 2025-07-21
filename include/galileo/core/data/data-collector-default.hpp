#ifndef __galileo_core_data_data_collector_default_hpp__
#define __galileo_core_data_data_collector_default_hpp__

#include "galileo/core/data/data-collector-actuation.hpp"
#include "galileo/core/data/data-collector-base.hpp"
#include "galileo/core/data/data-collector-contacts.hpp"
#include "galileo/core/data/data-collector-joint.hpp"
#include "galileo/core/data/data-collector-multibody.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct DataCollectorDefaultTpl : public DataCollectorBase<DataCollectorDefaultTpl<PhaseSpec, ContactCollectionTpl>>,
                                     public MultibodyDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                                     public ActuationDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                                     public JointDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                                     public ContactDataMixinTpl<DataCollectorDefaultTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec, ContactCollectionTpl>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using ContactManagerMeta_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using ContactModelManager_t = typename traits<ContactManagerMeta_t>::ModelManager_t;
        using ContactDataManager_t = typename traits<ContactManagerMeta_t>::DataManager_t;

        using RobotData_t = typename PS::RobotData_t;
        using ActuationData_t = typename PS::ActuationData_t;
        using JointData_t = JointDataTpl<PS>;
        using ContactData_t = ContactDataManager_t;

        DataCollectorDefaultTpl(RobotData_t *robot_,
                                ActuationData_t *actuation_,
                                JointData_t *joint_,
                                ContactData_t *contacts_) : MultibodyDataMixinTpl<DataCollectorDefaultTpl<PS, ContactCollectionTpl>, PS>(robot_),
                                                                ActuationDataMixinTpl<DataCollectorDefaultTpl<PS, ContactCollectionTpl>, PS>(actuation_),
                                                                JointDataMixinTpl<DataCollectorDefaultTpl<PS, ContactCollectionTpl>, PS>(joint_),
                                                                ContactDataMixinTpl<DataCollectorDefaultTpl<PS, ContactCollectionTpl>, PS, ContactCollectionTpl>(contacts_)

        {
        }

    }; // struct DataCollectorDefaultTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_default_hpp__
