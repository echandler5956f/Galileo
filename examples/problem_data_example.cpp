
#include <string>
#include <pinocchio/parsers/urdf.hpp>
#include <Eigen/Core>

using namespace galileo;

typedef double Scalar;

typedef pinocchio::ModelTpl<Scalar> Model;
typedef Model::Data Data;

typedef Model::ConfigVectorType ConfigVector;

const std::string huron_location = "resources/urdf/huron_cheat.urdf";
const int num_ees = 2;
const std::string end_effector_names[] = {"l_foot_v_ft_link", "r_foot_v_ft_link"};

void buildRobot(model::LeggedBody &bot)
{
    std::vector<std::string> ee_name_vect;
    ee_name_vect.resize(num_ees);
    for (int i = 0; i < num_ees; ++i)
    {
        ee_name_vect[i] = end_effector_names[i];
    }
    pinocchio::urdf::buildModel(huron_location, bot);

    bot.setEndEffectors(ee_name_vect);

    bot.GenerateContactCombination();
}

std::shared_ptr<environment::EnvironmentSurfaces> defineEnvironmentSurfaces()
{
    std::shared_ptr<environment::EnvironmentSurfaces> surfaces = std::make_shared<environment::EnvironmentSurfaces>();

    surfaces->push_back(environment::CreateInfiniteGround());

    return surfaces;
}

variables::States defineStates()
{
    variables::States states(19, 18);
    return states;
}

contact::ContactSequence defineContactSequence(variables::States states)
{
    contact::ContactSequence contact_sequence;

    contact::ContactMode initial_mode;
    initial_mode.combination_definition = bot.getContactCombination(0b11);
    initial_mode.contact_surfaces = {0, 0};

    contact::ContactMode intermediate_mode;
    intermediate_mode.combination_definition = bot.getContactCombination(0b10);
    intermediate_mode.contact_surfaces = {0, environment::NO_SURFACE};

    contact::ContactMode final_mode = initial_mode;

    std::cout << "Set the modes" << std::endl;

    // A contact sequence has timing and knot metadata
    contact::ContactSequence *contact_sequence = new contact::ContactSequence(num_ees);
    contact_sequence->addPhase(initial_mode, 100, 0.2);
    contact_sequence->addPhase(intermediate_mode, 100, 0.5);
    contact_sequence->addPhase(final_mode, 100, 0.3);

    return contact_sequence;
}

int main()
{
    model::LeggedBody bot;
    defineRobot(bot);

    std::shared_ptr<environment::EnvironmentSurfaces> surfaces = defineEnvironmentSurfaces();

    model::legged::LeggedProblemData;

    LeggedProblemData problem_data;

    return 0;
}