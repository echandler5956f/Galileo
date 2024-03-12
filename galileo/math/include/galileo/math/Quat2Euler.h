#pragma once

#include <Eigen/Dense>

namespace galileo
{
    namespace math
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

        void threeaxisrot(const Eigen::VectorXd &r11, const Eigen::VectorXd &r12, const Eigen::VectorXd &r21, const Eigen::VectorXd &r31, const Eigen::VectorXd &r32, Eigen::MatrixXd &res)
        {
            res.row(0) = r31.array().binaryExpr(r32.array(), [](double a, double b)
                                                { return atan2(a, b); });
            res.row(1) = r21.array().unaryExpr([](double a)
                                               { return asin(a); });
            res.row(2) = r11.array().binaryExpr(r12.array(), [](double a, double b)
                                                { return atan2(a, b); });
        }

        void twoaxisrot(const Eigen::VectorXd &r11, const Eigen::VectorXd &r12, const Eigen::VectorXd &r21, const Eigen::VectorXd &r31, const Eigen::VectorXd &r32, Eigen::MatrixXd &res)
        {
            res.row(0) = r11.array().binaryExpr(r12.array(), [](double a, double b)
                                                { return atan2(a, b); });
            res.row(1) = r21.array().unaryExpr([](double a)
                                               { return acos(a); });
            res.row(2) = r31.array().binaryExpr(r32.array(), [](double a, double b)
                                                { return atan2(a, b); });
        }

        Eigen::MatrixXd quaternion2Euler(const Eigen::MatrixXd &quat, RotSeq rotSeq)
        {
            Eigen::MatrixXd euler(3, quat.cols());

            Eigen::VectorXd r11, r12, r21, r31, r32;
            switch (rotSeq)
            {
            case zyx:
                r11 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(3).array() * quat.row(0).array());
                r12 = quat.row(3).array() * quat.row(3).array() - quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array();
                r21 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(3).array() * quat.row(1).array());
                r31 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(3).array() * quat.row(0).array());
                r32 = quat.row(3).array() * quat.row(3).array() + quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case zyz:
                r11 = 2 * (quat.row(2).array() * quat.row(3).array() - quat.row(0).array() * quat.row(1).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(1).array() * quat.row(3).array() + quat.row(0).array() * quat.row(2).array());
                r31 = 2 * (quat.row(2).array() * quat.row(3).array() - quat.row(0).array() * quat.row(1).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case zxy:
                r11 = -2 * (quat.row(1).array() * quat.row(3).array() - quat.row(0).array() * quat.row(2).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(2).array() * quat.row(3).array() + quat.row(0).array() * quat.row(1).array());
                r31 = -2 * (quat.row(1).array() * quat.row(3).array() - quat.row(0).array() * quat.row(2).array());
                r32 = quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case zxz:
                r11 = 2 * (quat.row(2).array() * quat.row(3).array() + quat.row(0).array() * quat.row(1).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                r21 = -2 * (quat.row(1).array() * quat.row(3).array() - quat.row(0).array() * quat.row(2).array());
                r31 = 2 * (quat.row(2).array() * quat.row(3).array() + quat.row(0).array() * quat.row(1).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case yxz:
                r11 = 2 * (quat.row(1).array() * quat.row(3).array() + quat.row(0).array() * quat.row(2).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                r21 = -2 * (quat.row(2).array() * quat.row(3).array() - quat.row(0).array() * quat.row(1).array());
                r31 = 2 * (quat.row(1).array() * quat.row(3).array() + quat.row(0).array() * quat.row(2).array());
                r32 = quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case yxy:
                r11 = 2 * (quat.row(1).array() * quat.row(2).array() - quat.row(0).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(0).array() * quat.row(3).array());
                r31 = 2 * (quat.row(0).array() * quat.row(1).array() - quat.row(2).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case yzx:
                r11 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(0).array() * quat.row(3).array());
                r31 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case yzy:
                r11 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(0).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r31 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(0).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case xyz:
                r11 = 2 * (quat.row(0).array() * quat.row(1).array() + quat.row(2).array() * quat.row(3).array());
                r12 = quat.row(3).array() * quat.row(3).array() - quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array();
                r21 = -2 * (quat.row(1).array() * quat.row(3).array() - quat.row(0).array() * quat.row(2).array());
                r31 = 2 * (quat.row(0).array() * quat.row(1).array() + quat.row(2).array() * quat.row(3).array());
                r32 = quat.row(3).array() * quat.row(3).array() + quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case xyx:
                r11 = 2 * (quat.row(0).array() * quat.row(1).array() - quat.row(2).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(0).array() * quat.row(1).array() + quat.row(2).array() * quat.row(3).array());
                r31 = 2 * (quat.row(1).array() * quat.row(2).array() - quat.row(0).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case xzy:
                r11 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(0).array() * quat.row(1).array() + quat.row(2).array() * quat.row(3).array());
                r31 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case xzx:
                r11 = 2 * (quat.row(0).array() * quat.row(2).array() + quat.row(1).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                r21 = -2 * (quat.row(0).array() * quat.row(1).array() - quat.row(2).array() * quat.row(3).array());
                r31 = 2 * (quat.row(0).array() * quat.row(2).array() + quat.row(1).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;
            default:
                printf("Rotation sequence not supported yet\n");
                break;
            }
            return euler;
        }
    }
}