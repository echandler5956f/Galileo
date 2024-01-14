#include "galileo/legged-model/LeggedBody.h"
#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/opt/PseudospectralSegment.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"

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

const int num_ees = 2;
const std::string end_effector_names[] = {"l_foot_v_ft_link", "r_foot_v_ft_link"};

/*So, for some reason, pinocchio's integrate function returns a slightly different SX matrix than than the python version.
When I actually evaluate the function though, it seems to return the correct result for both versions.
I have tried everything I can think of to figure out the discrepency, but I have not found the issue.
So, I just copied over the output from integrate for huron into this function here.
This is obviously a big problem for portability, but it at least lets us test if the rest of the framework is working.*/
SX custom_fint(SX x, SX dx, SX dt)
{
    SX x_12 = x(12);
    SX x_13 = x(13);
    SX x_14 = x(14);
    SX x_15 = x(15);
    SX x_16 = x(16);
    SX x_17 = x(17);
    SX x_18 = x(18);
    SX x_19 = x(19);
    SX x_20 = x(20);
    SX x_21 = x(21);
    SX x_22 = x(22);
    SX x_23 = x(23);
    SX x_24 = x(24);
    SX x_25 = x(25);
    SX x_26 = x(26);
    SX x_27 = x(27);
    SX x_28 = x(28);
    SX x_29 = x(29);
    SX x_30 = x(30);
    SX dx_12 = dx(12);
    SX dx_13 = dx(13);
    SX dx_14 = dx(14);
    SX dx_15 = dx(15);
    SX dx_16 = dx(16);
    SX dx_17 = dx(17);
    SX dx_18 = dx(18);
    SX dx_19 = dx(19);
    SX dx_20 = dx(20);
    SX dx_21 = dx(21);
    SX dx_22 = dx(22);
    SX dx_23 = dx(23);
    SX dx_24 = dx(24);
    SX dx_25 = dx(25);
    SX dx_26 = dx(26);
    SX dx_27 = dx(27);
    SX dx_28 = dx(28);
    SX dx_29 = dx(29);

    SX c1 = (dx_12 * dt);
    SX c2 = (dx_15 * dt);
    SX c3 = (dx_16 * dt);
    SX c4 = (dx_17 * dt);
    SX c5 = 4.93038e-32;
    SX c6 = ((sq(c2) + (sq(c3) + sq(c4))) + c5);
    SX c7 = sqrt(c6);
    SX c8 = 0.00012207;
    SX c9 = (c7 < c8);
    SX c10 = 0.5;
    SX c11 = 24;
    SX c12 = 1;
    SX c13 = (if_else(c9, (c10 - (c6 / c11)), 0) + if_else((!c9), ((c12 - cos(c7)) / c6), 0));
    SX c14 = (dx_14 * dt);
    SX c15 = (dx_13 * dt);
    SX c16 = (c7 < c8);
    SX c17 = 120;
    SX c18 = (if_else(c16, (0.166667 - (c6 / c17)), 0) + if_else((!c16), (((c7 - sin(c7)) / c6) / c7), 0));
    SX c19 = ((c2 * c15) - (c3 * c1));
    SX c20 = ((c4 * c1) - (c2 * c14));
    SX c21 = ((c1 + (c13 * ((c3 * c14) - (c4 * c15)))) + (c18 * ((c3 * c19) - (c4 * c20))));
    SX c22 = ((c3 * c14) - (c4 * c15));
    SX c23 = ((c14 + (c13 * ((c2 * c15) - (c3 * c1)))) + (c18 * ((c2 * c20) - (c3 * c22))));
    SX c24 = ((c15 + (c13 * ((c4 * c1) - (c2 * c14)))) + (c18 * ((c4 * c22) - (c2 * c19))));
    SX c25 = ((x_16 * c23) - (x_17 * c24));
    SX c26 = (c25 + c25);
    SX c27 = ((x_15 * c24) - (x_16 * c21));
    SX c28 = (c27 + c27);
    SX c29 = ((x_17 * c21) - (x_15 * c23));
    SX c30 = (c29 + c29);
    SX c31 = (sq(c2) + (sq(c3) + sq(c4)));
    SX c32 = (c8 < c31);
    SX c33 = sqrt((c31 + c5));
    SX c34 = (c10 * c33);
    SX c35 = sin(c34);
    SX c36 = (c31 / 4);
    SX c37 = (c10 * ((c12 - (c36 / 6)) + (sq(c36) / c17)));
    SX c38 = (if_else(c32, (c35 * (c2 / c33)), 0) + if_else((!c32), (c37 * c2), 0));
    SX c39 = (c8 < c31);
    SX c40 = 2;
    SX c41 = (if_else(c39, cos(c34), 0) + if_else((!c39), ((c12 - (c36 / c40)) + (sq(c36) / c11)), 0));
    SX c42 = (c8 < c31);
    SX c43 = (if_else(c42, (c35 * (c4 / c33)), 0) + if_else((!c42), (c37 * c4), 0));
    SX c44 = (c8 < c31);
    SX c45 = (if_else(c44, (c35 * (c3 / c33)), 0) + if_else((!c44), (c37 * c3), 0));
    SX c46 = ((((x_18 * c38) + (x_15 * c41)) + (x_16 * c43)) - (x_17 * c45));
    SX c47 = ((((x_18 * c45) + (x_16 * c41)) + (x_17 * c38)) - (x_15 * c43));
    SX c48 = ((((x_18 * c43) + (x_17 * c41)) + (x_15 * c45)) - (x_16 * c38));
    SX c49 = ((((x_18 * c41) - (x_15 * c38)) - (x_16 * c45)) - (x_17 * c43));
    SX c50 = (((c46 * x_15) + (c47 * x_16)) + ((c48 * x_17) + (c49 * x_18)));
    SX c51 = 0;
    SX c52 = (c50 < c51);
    SX c53 = (if_else(c52, (-c46), 0) + if_else((!c52), c46, 0));
    SX c54 = (c50 < c51);
    SX c55 = (if_else(c54, (-c47), 0) + if_else((!c54), c47, 0));
    SX c56 = (c50 < c51);
    SX c57 = (if_else(c56, (-c48), 0) + if_else((!c56), c48, 0));
    SX c58 = (c50 < c51);
    SX c59 = (if_else(c58, (-c49), 0) + if_else((!c58), c49, 0));
    SX c60 = ((3 - ((sq(c53) + sq(c55)) + (sq(c57) + sq(c59)))) / c40);

    return SX::vertcat(
        {(((c21 + (x_18 * c26)) + ((x_16 * c28) - (x_17 * c30))) + x_12),
         (((c24 + (x_18 * c30)) + ((x_17 * c26) - (x_15 * c28))) + x_13),
         (((c23 + (x_18 * c28)) + ((x_15 * c30) - (x_16 * c26))) + x_14),
         (c53 * c60),
         (c55 * c60),
         (c57 * c60),
         (c59 * c60),
         (x_19 + (dx_18 * dt)),
         (x_20 + (dx_19 * dt)),
         (x_21 + (dx_20 * dt)),
         (x_22 + (dx_21 * dt)),
         (x_23 + (dx_22 * dt)),
         (x_24 + (dx_23 * dt)),
         (x_25 + (dx_24 * dt)),
         (x_26 + (dx_25 * dt)),
         (x_27 + (dx_26 * dt)),
         (x_28 + (dx_27 * dt)),
         (x_29 + (dx_28 * dt)),
         (x_30 + (dx_29 * dt))});
}