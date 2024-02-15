#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/opt/TrajectoryOpt.h"

using namespace galileo;
using namespace legged;
using namespace opt;
using namespace constraints;

const std::string go1_location = "../resources/go1/urdf/go1_cheat.urdf";

const int num_ees = 4;
const std::string end_effector_names[] = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};