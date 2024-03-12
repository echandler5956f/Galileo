#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/tools/MeshcatInterface.h"
#include "galileo/tools/GNUPlotInterface.h"
#include "/home/quant/ros_ws/src/Galileo/galileo/legged-model/include/galileo/GalileoToOCS2Conversions.h"

using namespace galileo;
using namespace legged;
using namespace opt;
using namespace constraints;
using namespace tools;

const std::string go1_location = "../resources/go1/urdf/go1.urdf";

const int num_ees = 4;
const std::string end_effector_names[] = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};