#include "../model/include/LeggedBody.h"
#include "../model/include/ContactSequence.h"
#include "../model/include/EnvironmentSurfaces.h"

#include "../variables/include/States.h"
#include "../variables/include/TrajectoryGeneration.h"

#include <string>

#include <pinocchio/parsers/urdf.hpp>

using namespace acro;

const std::string huron_location = "resources/urdf/huron_cheat.urdf";
const int num_ees = 2;
const std::string end_effector_names[] = {"l_foot_v_ft_link", "r_foot_v_ft_link"};

void defineRobot(model::LeggedBody &bot)
{
    std::vector<std::string> ee_name_vect;
    ee_name_vect.resize(num_ees);
    for (int i = 0; i < num_ees; ++i)
    {
        ee_name_vect[i] = end_effector_names[i];
    }
    pinocchio::urdf::buildModel(huron_location, bot);
    std::cout << "Built model" << std::endl;
    bot.setEndEffectors(ee_name_vect);
    std::cout << "Set the end effectors" << std::endl;
    bot.GenerateContactCombination();
}

int main(int argc, char**argv)
{
    model::LeggedBody bot;
    defineRobot(bot);
    std::cout << "Set the robot" << std::endl;

    //Create environment Surfaces
    std::shared_ptr<environment::EnvironmentSurfaces> surfaces = std::make_shared<environment::EnvironmentSurfaces>();
    surfaces->push_back(environment::CreateInfiniteGround());
    std::cout << "Set the ground" << std::endl;

    // variables::Trajectory traj;
    // //  Initializes dynamics
    // traj.setModel(bot);
    
    
    // Target defines parameter for defining the cost
    // problemSetup defines parameters for setting up the problem, such as initial state
    variables::States* state_definition = new variables::States(19,18);
    Eigen::VectorXd target_state_vector;
    Eigen::VectorXd initial_state_vector;

    std::cout << "Set the states" << std::endl;

    contact::ContactMode initial_mode;
    initial_mode.combination_definition = bot.getContactCombination(0b11);
    initial_mode.contact_surfaces = {0, 0};

    
    contact::ContactMode intermediate_mode;
    intermediate_mode.combination_definition = bot.getContactCombination(0b10);
    intermediate_mode.contact_surfaces = {0,environment::NO_SURFACE};

    contact::ContactMode final_mode = initial_mode;
    
    std::cout << "Set the modes" << std::endl;
    
    // A contact sequence has timing and knot metadata
    contact::ContactSequence *contact_sequence = new contact::ContactSequence(num_ees);
    contact_sequence->addPhase(initial_mode, 100, 0.2);
    contact_sequence->addPhase(intermediate_mode, 100, 0.5);
    contact_sequence->addPhase(final_mode, 100, 0.3);


    std::cout << "Set the contact sequence" << std::endl;


    variables::Target<> target(target_state_vector, state_definition);
    variables::InitialCondition<> initial_condition(initial_state_vector, state_definition, initial_mode);
    variables::ProblemSetup<> problem_setup(initial_condition, contact_sequence);

    std::cout << "Finished setting up the problem" << std::endl;

    // //  Check the phases.
    // //  Check the contact combinations in these phases.
    contact::ContactSequence::Phase phase_0; 
    contact::ContactSequence::CONTACT_SEQUENCE_ERROR error_status_0;
    problem_setup.contact_sequence->getPhaseAtTime(0, phase_0, error_status_0);
    
    std::cout << "phase 0 time -- " << phase_0.time_value << std::endl;

    // std::cout << "phase 0 contacts -- ";
    // for(auto& end_effector_contact : phase_0.mode.combination_definition){
    //     if(end_effector_contact.second){
    //        std::cout << end_effector_contact.first << ",  ";
    //     }
    // }
    // std::cout << std::endl;

    // contact::ContactSequence::Phase phase_1; 
    // problem_setup.contact_sequence->getPhaseAtTime(phase_0.time_value + 0.01, phase_1, error_status_0);
        
    // std::cout << "phase 1 time -- " << phase_1.time_value << std::endl;

    // std::cout << "phase 1 contacts -- ";
    // for(auto& end_effector_contact : phase_1.mode.combination_definition){
    //     if(end_effector_contact.second){
    //        std::cout << end_effector_contact.first << ",  ";
    //     }
    // }
    // std::cout << std::endl;


    // contact::ContactSequence::Phase phase_2; 
    // problem_setup.contact_sequence->getPhaseAtTime(phase_0.time_value + phase_1.time_value + 0.01, phase_2, error_status_0);
        
    // std::cout << "phase 1 time -- " << phase_2.time_value << std::endl;

    // std::cout << "phase 1 contacts -- ";
    // for(auto& end_effector_contact : phase_2.mode.combination_definition){
    //     if(end_effector_contact.second){
    //        std::cout << end_effector_contact.first << ",  ";
    //     }
    // }
    // std::cout << std::endl;

    return 0;
}