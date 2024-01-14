#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/ContactSequence.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"

#include "galileo/legged-model/LeggedRobotStates.h"

#include <string>

#include <pinocchio/parsers/urdf.hpp>

using namespace galileo;
using namespace legged;

const std::string huron_location = "../resources/urdf/huron_cheat.urdf";
const int num_ees = 2;
const std::string end_effector_names[] = {"l_foot_v_ft_link", "r_foot_v_ft_link"};

int main()
{
    legged::LeggedBody<double> bot(huron_location, num_ees, end_effector_names);
    std::cout << "Set the robot" << std::endl;

    // Create environment Surfaces
    std::shared_ptr<environment::EnvironmentSurfaces> surfaces = std::make_shared<environment::EnvironmentSurfaces>();
    surfaces->push_back(environment::CreateInfiniteGround());
    std::cout << "Set the ground" << std::endl;

    // opt::Trajectory traj;
    // //  Initializes dynamics
    // traj.setModel(bot);

    // Target defines parameter for defining the cost
    // problemSetup defines parameters for setting up the problem, such as initial state
    std::shared_ptr<opt::LeggedRobotStates> si = std::make_shared<opt::LeggedRobotStates>(std::vector<int>{bot.model.nq, bot.model.nv});
    Eigen::VectorXd target_state_vector;
    Eigen::VectorXd initial_state_vector;

    std::cout << "Set the states" << std::endl;

    contact::ContactMode initial_mode;
    initial_mode.combination_definition = bot.getContactCombination(0b11);
    initial_mode.contact_surfaces = {0, 0};

    contact::ContactMode intermediate_mode;
    auto combination = bot.getContactCombination(0b10);
    intermediate_mode.combination_definition = combination;
    intermediate_mode.contact_surfaces.resize(combination.size());

    // Set the surfaces for the end effectors that are in contact
    int i = 0;
    for (auto combo : combination)
    {
        if (combo.second) // If the end effector is in contact
        {
            // Set the surface for this end effector
            intermediate_mode.contact_surfaces[i] = 0;
        }
        else
        {
            intermediate_mode.contact_surfaces[i] = environment::NO_SURFACE;
        }
    }

    contact::ContactMode final_mode = initial_mode;

    std::cout << "Set the modes" << std::endl;

    // A contact sequence has timing and knot metadata
    contact::ContactSequence *contact_sequence = new contact::ContactSequence(num_ees);
    contact_sequence->addPhase(initial_mode, 100, 0.2);
    contact_sequence->addPhase(intermediate_mode, 100, 0.5);
    contact_sequence->addPhase(final_mode, 100, 0.3);

    std::cout << "Set the contact sequence" << std::endl;

    // opt::Target<> target(target_state_vector, state_definition);
    // opt::InitialCondition<> initial_condition(initial_state_vector, state_definition, initial_mode);
    // opt::ProblemSetup<> problem_setup(initial_condition, contact_sequence);

    std::cout << "Finished setting up the problem" << std::endl;

    // //  Check the phases.
    // //  Check the contact combinations in these phases.
    contact::ContactSequence::Phase phase_0;
    contact::ContactSequence::CONTACT_SEQUENCE_ERROR error_status_0;
    contact_sequence->getPhaseAtTime(0, phase_0, error_status_0);

    std::cout << "phase 0 time -- " << phase_0.time_value << std::endl;

    std::cout << "phase 0 contacts -- ";
    for (auto &end_effector_contact : phase_0.mode.combination_definition)
    {
        if (end_effector_contact.second)
        {
            std::cout << end_effector_contact.first << ",  ";
        }
    }
    std::cout << std::endl;

    contact::ContactSequence::Phase phase_1;
    contact_sequence->getPhaseAtTime(phase_0.time_value + 0.01, phase_1, error_status_0);

    std::cout << "phase 1 time -- " << phase_1.time_value << std::endl;

    std::cout << "phase 1 contacts -- ";
    for (auto &end_effector_contact : phase_1.mode.combination_definition)
    {
        if (end_effector_contact.second)
        {
            std::cout << end_effector_contact.first << ",  ";
        }
    }
    std::cout << std::endl;

    contact::ContactSequence::Phase phase_2;
    contact_sequence->getPhaseAtTime(phase_0.time_value + phase_1.time_value + 0.01, phase_2, error_status_0);

    std::cout << "phase 2 time -- " << phase_2.time_value << std::endl;

    std::cout << "phase 2 contacts -- ";
    for (auto &end_effector_contact : phase_2.mode.combination_definition)
    {
        if (end_effector_contact.second)
        {
            std::cout << end_effector_contact.first << ",  ";
        }
    }
    std::cout << std::endl;

    return 0;
}