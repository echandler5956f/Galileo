#include "galileo/model/LeggedBody.h"
#include "galileo/opt/TrajectoryOpt.h"
#include <string>

#include <pinocchio/parsers/urdf.hpp>
#include <Eigen/Dense>

#include <pinocchio/autodiff/casadi.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>

#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/centroidal-derivatives.hpp>

using namespace galileo;
using namespace casadi;
using namespace std;

typedef double Scalar;
typedef SX ADScalar;

typedef pinocchio::ModelTpl<Scalar> Model;
typedef Model::Data Data;

typedef pinocchio::ModelTpl<ADScalar> ADModel;
typedef ADModel::Data ADData;

typedef Model::ConfigVectorType ConfigVector;
typedef Model::TangentVectorType TangentVector;

typedef ADModel::ConfigVectorType ConfigVectorAD;
typedef ADModel::TangentVectorType TangentVectorAD;

const string huron_location = "../resources/urdf/huron_cheat.urdf";