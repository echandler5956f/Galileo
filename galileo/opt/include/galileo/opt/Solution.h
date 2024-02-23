#pragma once

#include <Eigen/Dense>
#include <string>

namespace galileo
{
    namespace opt
    {
        using tuple_size_t = std::tuple<size_t, size_t>;

        struct constraint_metadata_t
        {
            /**
             * @brief Name of the constraint.
             *
             */
            std::string name;

            /**
             * @brief A vector of index ranges which group certain constraint rows together to make plotting easier. An empty vector means to include all constraint rows in the same plot.
             *
             */
            std::vector<tuple_size_t> plot_groupings;

            /**
             * @brief Titles of the plots.
             *
             */
            std::vector<std::string> plot_titles;

            /**
             * @brief Names of the plots on the legend.
             *
             */
            std::vector<std::vector<std::string>> plot_names;
        };

        struct constraint_evaluations_t
        {
            Eigen::VectorXd times;
            Eigen::MatrixXd evaluation;
            Eigen::MatrixXd lower_bounds;
            Eigen::MatrixXd upper_bounds;
            constraint_metadata_t metadata;
        };

        struct solution_t
        {
            solution_t() {}
            solution_t(Eigen::VectorXd times_) { this->times = times_; }
            Eigen::VectorXd times;
            Eigen::MatrixXd state_result;
            Eigen::MatrixXd input_result;
        };
    }
}