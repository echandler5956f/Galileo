#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/legged-model/LeggedInterface.h"
#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/tools/MeshcatInterface.h"
#include "galileo/tools/GNUPlotInterface.h"
#include "galileo/math/Quat2Euler.h"

using namespace galileo;
using namespace legged;
using namespace opt;
using namespace constraints;
using namespace tools;
using namespace math;

const std::string go1_location = "../resources/go1/urdf/go1.urdf";

const int num_ees = 4;
const std::vector<std::string> end_effector_names = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};