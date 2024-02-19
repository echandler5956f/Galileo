#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/opt/TrajectoryOpt.h"

using namespace galileo;
using namespace legged;
using namespace opt;
using namespace constraints;

const std::string atlas_location = "../resources/atlas/atlas_minimal_contact.urdf";

const int num_ees = 2;
const std::string end_effector_names[] = {"l_foot", "r_foot"};