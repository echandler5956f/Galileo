#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/legged-model/LeggedInterface.h"
#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/tools/MeshcatInterface.h"
#include "galileo/tools/GNUPlotInterface.h"
#include "galileo/math/Quat2Euler.h"

#include "galileo/legged-model/LeggedModelHelpers.h"

using namespace galileo;
using namespace legged;
using namespace opt;
using namespace constraints;
using namespace tools;
using namespace math;

const std::string robot_location = "../resources/atlas/urdf/atlas.urdf";
const std::string solver_parameter_location = "../resources/atlas/Parameters/solver_parameters.txt";
const std::string problem_parameter_location = "../resources/atlas/Parameters/problem_parameters.txt";