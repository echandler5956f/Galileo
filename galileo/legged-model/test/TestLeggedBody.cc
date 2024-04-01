#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>
#include "galileo/legged-model/LeggedBody.h"

using namespace galileo::legged;

TEST(LeggedBodyTest, LeggedBodyConstruction)
{
    std::string location = "../resources/go1/urdf/go1.urdf";;
    std::vector<std::string> end_effector_names = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};
    LeggedBody robot(location, end_effector_names);

    int nq = 19;
    int nv = 18;
    EXPECT_EQ(robot.si->nq, nq);
    EXPECT_EQ(robot.si->nv, nv);
    EXPECT_EQ(robot.si->nx, robot.si->nh + robot.si->nq);
    EXPECT_EQ(robot.si->ndx, robot.si->ndh + robot.si->nv);
    EXPECT_EQ(robot.si->nvju, robot.si->nv - robot.si->nvb);
    EXPECT_EQ(robot.si->nu, robot.si->nF + robot.si->nvju);
    EXPECT_EQ(robot.si->nu_general, 6 + robot.si->nvju);
    EXPECT_EQ(robot.si->h_index, 0);
    EXPECT_EQ(robot.si->q_index, robot.si->h_index + robot.si->nh);
    EXPECT_EQ(robot.si->qj_index, robot.si->q_index + robot.si->nqb);
    EXPECT_EQ(robot.si->general_force_index, 0);
    EXPECT_EQ(robot.si->general_torque_index, robot.si->general_force_index + 3);
    EXPECT_EQ(robot.si->general_vju_index, robot.si->general_torque_index + 3);
}