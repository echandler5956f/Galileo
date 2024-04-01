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

const std::string atlas_location = "../resources/atlas/urdf/atlas.urdf";
const std::string atlas_parameter_location = "../resources/atlas/SolverParameters/solver_parameters.txt";
const std::vector<std::string> end_effector_names = {"l_foot", "r_foot"};