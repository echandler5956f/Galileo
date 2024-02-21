#pragma once

#include <Eigen/Dense>

namespace galileo
{
    namespace opt
    {
        struct constraint_evaluations_t
        {
            std::string name;
            Eigen::Matrix3Xd evaluation_and_bounds;
        };

        struct solution_t
        {
            solution_t(Eigen::VectorXd times_) { this->times = times_; }
            Eigen::VectorXd times;
            Eigen::MatrixXd state_result;
            Eigen::MatrixXd input_result;
        };
    }
}