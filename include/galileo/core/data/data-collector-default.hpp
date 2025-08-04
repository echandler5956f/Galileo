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
                                     public JointDataMixinTpl<DataCollectorContactTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>
    {
        using PS = PhaseSpec;

        using RobotData_t = typename PS::RobotData_t;
        using ActuationData_t = typename PS::ActuationData_t;
        using JointData_t = JointDataTpl<PS>;

        DataCollectorContactTpl(std::shared_ptr<RobotData_t> robot_,
                                std::shared_ptr<ActuationData_t> actuation_,
                                std::shared_ptr<JointData_t> joint_) : MultibodyDataMixinTpl<DataCollectorContactTpl<PS, ContactCollectionTpl>, PS>(robot_),
                                                                       ActuationDataMixinTpl<DataCollectorContactTpl<PS, ContactCollectionTpl>, PS>(actuation_),
                                                                       JointDataMixinTpl<DataCollectorContactTpl<PS, ContactCollectionTpl>, PS>(joint_)
        {
        }

    }; // struct DataCollectorContactTpl

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct DataCollectorImpulseTpl : public DataCollectorBase<DataCollectorImpulseTpl<PhaseSpec, ImpulseCollectionTpl>>,
                                     public MultibodyDataMixinTpl<DataCollectorImpulseTpl<PhaseSpec, ImpulseCollectionTpl>, PhaseSpec>
    {
        using PS = PhaseSpec;

        using ImpulseManagerMeta_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using ImpulseModelManager_t = typename traits<ImpulseManagerMeta_t>::ModelManager_t;
        using ImpulseDataManager_t = typename traits<ImpulseManagerMeta_t>::DataManager_t;

        using RobotData_t = typename PS::RobotData_t;
        using ImpulseData_t = ImpulseDataManager_t;

        DataCollectorImpulseTpl(std::shared_ptr<RobotData_t> robot_) : MultibodyDataMixinTpl<DataCollectorImpulseTpl<PS, ImpulseCollectionTpl>, PS>(robot_)
        {
        }

    }; // struct DataCollectorImpulseTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_default_hpp__
