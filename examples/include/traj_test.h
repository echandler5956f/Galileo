#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/opt/TrajectoryOpt.h"

using namespace galileo;
using namespace legged;
using namespace opt;
using namespace constraints;

const std::string huron_location = "../resources/urdf/huron_cheat.urdf";

const int num_ees = 2;
const std::string end_effector_names[] = {"l_foot_v_ft_link", "r_foot_v_ft_link"};