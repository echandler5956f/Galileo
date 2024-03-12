#pragma once

#include <Eigen/Dense>

namespace galileo
{
    namespace legged
    {
        /*Taken from https://stackoverflow.com/questions/1031005/is-there-an-algorithm-for-converting-quaternion-rotations-to-euler-angle-rotatio*/
        enum RotSeq
        {
            zyx,
            zyz,
            zxy,
            zxz,
            yxz,
            yxy,
            yzx,
            yzy,
            xyz,
            xyx,
            xzy,
            xzx
        };

        /*TODO: Make atan2*/
void threeaxisrot(const Eigen::VectorXd &r11, const Eigen::VectorXd &r12, const Eigen::VectorXd &r21, const Eigen::VectorXd &r31, const Eigen::VectorXd &r32, Eigen::MatrixXd &res)
{
    res.row(0) = r31.array().binaryExpr(r32.array(), [](double a, double b) { return atan2(a, b); });
    res.row(1) = r21.array().unaryExpr([](double a) { return asin(a); });
    res.row(2) = r11.array().binaryExpr(r12.array(), [](double a, double b) { return atan2(a, b); });
}

void twoaxisrot(const Eigen::VectorXd &r11, const Eigen::VectorXd &r12, const Eigen::VectorXd &r21, const Eigen::VectorXd &r31, const Eigen::VectorXd &r32, Eigen::MatrixXd &res)
{
    res.row(0) = r11.array().binaryExpr(r12.array(), [](double a, double b) { return atan2(a, b); });
    res.row(1) = r21.array().unaryExpr([](double a) { return acos(a); });
    res.row(2) = r31.array().binaryExpr(r32.array(), [](double a, double b) { return atan2(a, b); });
}

        Eigen::MatrixXd quaternion2Euler(const Eigen::MatrixXd &quat, RotSeq rotSeq)
        {
            Eigen::MatrixXd euler(3, quat.cols());

            // switch (rotSeq)
            // {
            // case zyx:
                Eigen::VectorXd r11 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(3).array() * quat.row(0).array());
                Eigen::VectorXd r12 = quat.row(3).array() * quat.row(3).array() - quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array();
                Eigen::VectorXd r21 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(3).array() * quat.row(1).array());
                Eigen::VectorXd r31 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(3).array() * quat.row(0).array());
                Eigen::VectorXd r32 = quat.row(3).array() * quat.row(3).array() + quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                // break;
            // default:
            //     std::cout << "Rotation sequence not supported yet" << std::endl;
            //     break;
            // }
            return euler;
        }

        Eigen::MatrixXd quatToEulerZYX(const Eigen::MatrixXd &quat)
        {
            return quaternion2Euler(quat, zyx);
        }
    }
}