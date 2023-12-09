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
    std::cout << "Starting" << std::endl;
    defineRobot(bot);
    std::cout << "Finished" << std::endl;

    // //Create environment Surfaces
    // std::shared_ptr<environment::EnvironmentSurfaces> surfaces;
    // surfaces->push_back(environment::CreateInfiniteGround());

    // std::cout << "surfaces->push_back" << std::endl;

    // // variables::Trajectory traj;
    // // //  Initializes dynamics
    // // traj.setModel(bot);
    
    
    // // Target defines parameter for defining the cost
    // // problemSetup defines parameters for setting up the problem, such as initial state
    // variables::States state_definition;
    // Eigen::VectorXd target_state_vector;
    // Eigen::VectorXd initial_state_vector;

    // contact::ContactMode initial_mode;
    // initial_mode.combination_definition = bot.getContactCombination(0b11);
    // initial_mode.contact_surfaces = {0, 0};
    // std::cout << "Here4" << std::endl;

    
    // contact::ContactMode intermediate_mode;
    // intermediate_mode.combination_definition = bot.getContactCombination(0b10);
    // intermediate_mode.contact_surfaces = {0,environment::NO_SURFACE};

    // contact::ContactMode final_mode = initial_mode;

    // std::cout << "Here5" << std::endl;

    
    // // A contact sequence has timing and knot metadata
    // contact::ContactSequence *contact_sequence = new contact::ContactSequence(num_ees);
    // contact_sequence->addPhase(initial_mode, 100, 0.2);
    // contact_sequence->addPhase(intermediate_mode, 100, 0.5);
    // contact_sequence->addPhase(final_mode, 100, 0.3);

    // std::cout << "Here6" << std::endl;


    // variables::Target<> target(target_state_vector, state_definition);
    // variables::InitialCondition<> initial_condition(initial_state_vector, state_definition, initial_mode);
    // // variables::ProblemSetup<> problem_setup(initial_condition, contact_sequence);

    // std::cout << "Here7" << std::endl;
    return -1;
}