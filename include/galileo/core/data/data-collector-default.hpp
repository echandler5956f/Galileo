#ifndef __galileo_core_data_data_collector_default_hpp__
#define __galileo_core_data_data_collector_default_hpp__

#include "galileo/core/data/data-collector-actuation.hpp"
#include "galileo/core/data/data-collector-base.hpp"
#include "galileo/core/data/data-collector-contacts.hpp"
#include "galileo/core/data/data-collector-impulses.hpp"
#include "galileo/core/data/data-collector-joint.hpp"
#include "galileo/core/data/data-collector-multibody.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct DataCollectorContactTpl : public DataCollectorBase<DataCollectorContactTpl<PhaseSpec, ContactCollectionTpl>>,
                                     public MultibodyDataMixinTpl<DataCollectorContactTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                                     public ActuationDataMixinTpl<DataCollectorContactTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                                     public JointDataMixinTpl<DataCollectorContactTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                                     public ContactDataMixinTpl<DataCollectorContactTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec, ContactCollectionTpl>
    {
        using PS = PhaseSpec;

        using ContactManagerMeta_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using ContactModelManager_t = typename traits<ContactManagerMeta_t>::ModelManager_t;
        using ContactDataManager_t = typename traits<ContactManagerMeta_t>::DataManager_t;

        using RobotData_t = typename PS::RobotData_t;
        using ActuationData_t = typename PS::ActuationData_t;
        using JointData_t = JointDataTpl<PS>;
        using ContactData_t = ContactDataManager_t;

        DataCollectorContactTpl(RobotData_t *robot_,
                                ActuationData_t *actuation_,
                                JointData_t *joint_,
                                ContactData_t *contacts_) : MultibodyDataMixinTpl<DataCollectorContactTpl<PS, ContactCollectionTpl>, PS>(robot_),
                                                            ActuationDataMixinTpl<DataCollectorContactTpl<PS, ContactCollectionTpl>, PS>(actuation_),
                                                            JointDataMixinTpl<DataCollectorContactTpl<PS, ContactCollectionTpl>, PS>(joint_),
                                                            ContactDataMixinTpl<DataCollectorContactTpl<PS, ContactCollectionTpl>, PS, ContactCollectionTpl>(contacts_)

        {
        }

    }; // struct DataCollectorContactTpl

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct DataCollectorImpulseTpl : public DataCollectorBase<DataCollectorImpulseTpl<PhaseSpec, ImpulseCollectionTpl>>,
                                     public MultibodyDataMixinTpl<DataCollectorImpulseTpl<PhaseSpec, ImpulseCollectionTpl>, PhaseSpec>,
                                     public ImpulseDataMixinTpl<DataCollectorImpulseTpl<PhaseSpec, ImpulseCollectionTpl>, PhaseSpec, ImpulseCollectionTpl>
    {
        using PS = PhaseSpec;

        using ImpulseManagerMeta_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using ImpulseModelManager_t = typename traits<ImpulseManagerMeta_t>::ModelManager_t;
        using ImpulseDataManager_t = typename traits<ImpulseManagerMeta_t>::DataManager_t;

        using RobotData_t = typename PS::RobotData_t;
        using ImpulseData_t = ImpulseDataManager_t;

        DataCollectorImpulseTpl(RobotData_t *robot_,
                                ImpulseData_t *impulses_) : MultibodyDataMixinTpl<DataCollectorImpulseTpl<PS, ImpulseCollectionTpl>, PS>(robot_),
                                                            ImpulseDataMixinTpl<DataCollectorImpulseTpl<PS, ImpulseCollectionTpl>, PS, ImpulseCollectionTpl>(impulses_)

        {
        }

    }; // struct DataCollectorImpulseTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_default_hpp__
