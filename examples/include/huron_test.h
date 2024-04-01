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

const std::string huron_location = "../resources/huron/urdf/huron.urdf";
const std::string huron_parameter_location = "../resources/huron/SolverParameters/solver_parameters.txt";
const std::vector<std::string> end_effector_names = {"l_foot_v_ft_link", "r_foot_v_ft_link"};