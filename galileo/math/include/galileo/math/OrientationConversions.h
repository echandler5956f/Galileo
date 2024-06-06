#pragma once

#include "galileo/math/OrientationDefinition.h"
#include <Eigen/Dense>

namespace galileo
{
    namespace math
    {
        /*Modified from https://stackoverflow.com/questions/1031005/is-there-an-algorithm-for-converting-quaternion-rotations-to-euler-angle-rotatio*/

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

        Eigen::MatrixXd quaternion2Euler(const Eigen::MatrixXd &quat, OrientationDefinition output_euler_def)
        {
            Eigen::MatrixXd euler(3, quat.cols());

            Eigen::VectorXd r11, r12, r21, r31, r32;
            switch (output_euler_def)
            {
            case OrientationDefinition::EulerZYX:
                r11 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(3).array() * quat.row(0).array());
                r12 = quat.row(3).array() * quat.row(3).array() - quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array();
                r21 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(3).array() * quat.row(1).array());
                r31 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(3).array() * quat.row(0).array());
                r32 = quat.row(3).array() * quat.row(3).array() + quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerZYZ:
                r11 = 2 * (quat.row(2).array() * quat.row(3).array() - quat.row(0).array() * quat.row(1).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(1).array() * quat.row(3).array() + quat.row(0).array() * quat.row(2).array());
                r31 = 2 * (quat.row(2).array() * quat.row(3).array() - quat.row(0).array() * quat.row(1).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerZXY:
                r11 = -2 * (quat.row(1).array() * quat.row(3).array() - quat.row(0).array() * quat.row(2).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(2).array() * quat.row(3).array() + quat.row(0).array() * quat.row(1).array());
                r31 = -2 * (quat.row(1).array() * quat.row(3).array() - quat.row(0).array() * quat.row(2).array());
                r32 = quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerZXZ:
                r11 = 2 * (quat.row(2).array() * quat.row(3).array() + quat.row(0).array() * quat.row(1).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                r21 = -2 * (quat.row(1).array() * quat.row(3).array() - quat.row(0).array() * quat.row(2).array());
                r31 = 2 * (quat.row(2).array() * quat.row(3).array() + quat.row(0).array() * quat.row(1).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerYXZ:
                r11 = 2 * (quat.row(1).array() * quat.row(3).array() + quat.row(0).array() * quat.row(2).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                r21 = -2 * (quat.row(2).array() * quat.row(3).array() - quat.row(0).array() * quat.row(1).array());
                r31 = 2 * (quat.row(1).array() * quat.row(3).array() + quat.row(0).array() * quat.row(2).array());
                r32 = quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerYXY:
                r11 = 2 * (quat.row(1).array() * quat.row(2).array() - quat.row(0).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(0).array() * quat.row(3).array());
                r31 = 2 * (quat.row(0).array() * quat.row(1).array() - quat.row(2).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerYZX:
                r11 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(0).array() * quat.row(3).array());
                r31 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerYZY:
                r11 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(0).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r31 = 2 * (quat.row(1).array() * quat.row(2).array() + quat.row(0).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerXYZ:
                r11 = 2 * (quat.row(0).array() * quat.row(1).array() + quat.row(2).array() * quat.row(3).array());
                r12 = quat.row(3).array() * quat.row(3).array() - quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array();
                r21 = -2 * (quat.row(1).array() * quat.row(3).array() - quat.row(0).array() * quat.row(2).array());
                r31 = 2 * (quat.row(0).array() * quat.row(1).array() + quat.row(2).array() * quat.row(3).array());
                r32 = quat.row(3).array() * quat.row(3).array() + quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerXYX:
                r11 = 2 * (quat.row(0).array() * quat.row(1).array() - quat.row(2).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(0).array() * quat.row(1).array() + quat.row(2).array() * quat.row(3).array());
                r31 = 2 * (quat.row(1).array() * quat.row(2).array() - quat.row(0).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() + quat.row(3).array() * quat.row(3).array();
                twoaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerXZY:
                r11 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r12 = quat.row(0).array() * quat.row(0).array() - quat.row(1).array() * quat.row(1).array() + quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                r21 = 2 * (quat.row(0).array() * quat.row(1).array() + quat.row(2).array() * quat.row(3).array());
                r31 = -2 * (quat.row(0).array() * quat.row(2).array() - quat.row(1).array() * quat.row(3).array());
                r32 = quat.row(0).array() * quat.row(0).array() + quat.row(1).array() * quat.row(1).array() - quat.row(2).array() * quat.row(2).array() - quat.row(3).array() * quat.row(3).array();
                threeaxisrot(r11, r12, r21, r31, r32, euler);
                break;

            case OrientationDefinition::EulerXZX:
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

        Eigen::MatrixXd euler2Quat(const Eigen::MatrixXd &euler, OrientationDefinition input_euler_def)
        {
            Eigen::MatrixXd quat(4, euler.cols());
            Eigen::MatrixXd halved_angles = euler.array() / 2.0;

            Eigen::VectorXd s1 = halved_angles.row(0).array().sin();
            Eigen::VectorXd s2 = halved_angles.row(1).array().sin();
            Eigen::VectorXd s3 = halved_angles.row(2).array().sin();
            Eigen::VectorXd c1 = halved_angles.row(0).array().cos();
            Eigen::VectorXd c2 = halved_angles.row(1).array().cos();
            Eigen::VectorXd c3 = halved_angles.row(2).array().cos();

            switch (input_euler_def)
            {
            case OrientationDefinition::EulerZYX:
                quat.row(0) = (s1.array() * c2.array() * c3.array()) - (c1.array() * s2.array() * s3.array());
                quat.row(1) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(2) = (c1.array() * c2.array() * s3.array()) - (s1.array() * s2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) + (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerZYZ:
                quat.row(0) = (s1.array() * c2.array() * c3.array()) + (c1.array() * s2.array() * s3.array());
                quat.row(1) = (s1.array() * s2.array() * c3.array()) - (c1.array() * c2.array() * s3.array());
                quat.row(2) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) - (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerZXZ:
                quat.row(0) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(1) = (c1.array() * c2.array() * s3.array()) - (s1.array() * s2.array() * c3.array());
                quat.row(2) = (c1.array() * s2.array() * s3.array()) + (s1.array() * c2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) - (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerYXZ:
                quat.row(0) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(1) = (c1.array() * c2.array() * s3.array()) + (s1.array() * s2.array() * c3.array());
                quat.row(2) = (c1.array() * s2.array() * s3.array()) - (s1.array() * c2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) - (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerYXY:
                quat.row(0) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(1) = (c1.array() * c2.array() * s3.array()) + (s1.array() * s2.array() * c3.array());
                quat.row(2) = (c1.array() * s2.array() * s3.array()) - (s1.array() * c2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) - (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerYZX:
                quat.row(0) = (s1.array() * c2.array() * c3.array()) - (c1.array() * s2.array() * s3.array());
                quat.row(1) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(2) = (c1.array() * c2.array() * s3.array()) - (s1.array() * s2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) + (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerYZY:
                quat.row(0) = (s1.array() * c2.array() * c3.array()) + (c1.array() * s2.array() * s3.array());
                quat.row(1) = (s1.array() * s2.array() * c3.array()) - (c1.array() * c2.array() * s3.array());
                quat.row(2) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) - (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerXYZ:
                quat.row(0) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(1) = (c1.array() * c2.array() * s3.array()) + (s1.array() * s2.array() * c3.array());
                quat.row(2) = (c1.array() * s2.array() * s3.array()) - (s1.array() * c2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) - (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerXYX:
                quat.row(0) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(1) = (c1.array() * c2.array() * s3.array()) - (s1.array() * s2.array() * c3.array());
                quat.row(2) = (c1.array() * s2.array() * s3.array()) + (s1.array() * c2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) - (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerXZY:
                quat.row(0) = (s1.array() * c2.array() * c3.array()) - (c1.array() * s2.array() * s3.array());
                quat.row(1) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(2) = (c1.array() * c2.array() * s3.array()) - (s1.array() * s2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) + (s1.array() * s2.array() * s3.array());
                break;
            case OrientationDefinition::EulerXZX:
                quat.row(0) = (c1.array() * s2.array() * c3.array()) + (s1.array() * c2.array() * s3.array());
                quat.row(1) = (c1.array() * c2.array() * s3.array()) - (s1.array() * s2.array() * c3.array());
                quat.row(2) = (c1.array() * s2.array() * s3.array()) + (s1.array() * c2.array() * c3.array());
                quat.row(3) = (c1.array() * c2.array() * c3.array()) - (s1.array() * s2.array() * s3.array());
                break;
            default:
                printf("Rotation sequence not supported yet\n");
                break;
            }
            return quat;
        }

        Eigen::MatrixXd OrientationConversion(const Eigen::MatrixXd &orientation_input, OrientationDefinition orientation_def_input, OrientationDefinition orientation_def_output)
        {
            if (orientation_def_input == orientation_def_output)
            {
                return orientation_input;
            }
            else if (orientation_def_input == OrientationDefinition::Quaternion)
            {
                return quaternion2Euler(orientation_input, orientation_def_output);
            }
            else if (orientation_def_output == OrientationDefinition::Quaternion)
            {
                return euler2Quat(orientation_input, orientation_def_input);
            }
            else
            {
                printf("Conversion not supported yet\n");
                return Eigen::MatrixXd::Zero(0, 0);
            }
        }
    }
}