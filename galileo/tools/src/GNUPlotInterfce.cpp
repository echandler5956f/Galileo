#include "galileo/tools/GNUPlotInterface.h"

namespace galileo
{
    namespace tools
    {
        GNUPlotInterface::GNUPlotInterface(opt::solution_t solution_, std::vector<std::vector<opt::constraint_evaluations_t>> constraint_evaluations_)
        {
            this->solution.state_result = solution_.state_result.transpose();
            this->solution.input_result = solution_.input_result.transpose();
            this->solution.times = solution_.times;
            this->constraint_evaluations = constraint_evaluations_;
        }

        void GNUPlotInterface::PlotSolution(std::vector<tuple_size_t> state_groups, std::vector<tuple_size_t> input_groups,
                                            std::vector<std::string> state_title_names, std::vector<std::vector<std::string>> state_names,
                                            std::vector<std::string> input_title_names, std::vector<std::vector<std::string>> input_names)
        {
            std::vector<double> std_times(solution.times.data(), solution.times.data() + solution.times.size());

            for (size_t i = 0; i < state_groups.size(); ++i)
            {
                gp << "set term wxt " << window_index << "\n"; // Create a new window for each plot
                gp << "set xlabel 'Time'\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + state_title_names[i] + "'\n";
                std::string ss = "plot ";
                for (size_t j = 0; j < state_names[i].size(); ++j)
                {
                    if (j != 0)
                    {
                        ss += ", ";
                    }
                    ss += "'-' with linespoints linestyle " + std::to_string(j + 1) + " title '" + state_names[i][j] + "'";
                }
                ss += "\n";
                gp << ss;
                Eigen::MatrixXd block = solution.state_result.block(0, std::get<0>(state_groups[i]), solution.state_result.rows(), std::get<1>(state_groups[i]) - std::get<0>(state_groups[i]));
                for (Eigen::Index j = 0; j < block.cols(); ++j)
                {
                    Eigen::VectorXd colMatrix = block.col(j).matrix();
                    std::vector<double> std_col_vector(colMatrix.data(), colMatrix.data() + colMatrix.size());
                    gp.send1d(std::make_tuple(std_times, std_col_vector));
                }
                ++window_index;
            }

            for (size_t i = 0; i < input_groups.size(); ++i)
            {
                gp << "set term wxt " << window_index << "\n"; // Create a new window for each plot
                gp << "set xlabel 'Time'\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + input_title_names[i] + "'\n";
                std::string ss = "plot ";
                for (size_t j = 0; j < input_names[i].size(); ++j)
                {
                    if (j != 0)
                    {
                        ss += ", ";
                    }
                    ss += "'-' with linespoints linestyle " + std::to_string(j + 1) + " title '" + input_names[i][j] + "'";
                }
                ss += "\n";
                gp << ss;
                Eigen::MatrixXd block = solution.input_result.block(0, std::get<0>(input_groups[i]), solution.input_result.rows(), std::get<1>(input_groups[i]) - std::get<0>(input_groups[i]));
                for (Eigen::Index j = 0; j < block.cols(); ++j)
                {
                    Eigen::VectorXd colMatrix = block.col(j).matrix();
                    std::vector<double> std_col_vector(colMatrix.data(), colMatrix.data() + colMatrix.size());
                    gp.send1d(std::make_tuple(std_times, std_col_vector));
                }
                ++window_index;
            }
        }

        void GNUPlotInterface::PlotConstraints()
        {
            std::map<std::string, opt::constraint_evaluations_t> constraint_map;

            for (auto &cons_in_phase : constraint_evaluations)
            {
                for (auto &constraint : cons_in_phase)
                {
                    for (size_t i = 0; i < constraint.metadata.plot_titles.size(); ++i)
                    {
                        std::string key = constraint.metadata.name + constraint.metadata.plot_titles[i];

                        if (constraint_map.count(key) > 0)
                        {
                            // If the constraint already exists, concatenate the new data
                            constraint_map[key].evaluation.conservativeResize(constraint_map[key].evaluation.rows(), constraint_map[key].evaluation.cols() + constraint.evaluation.cols());
                            constraint_map[key].evaluation << constraint.evaluation;

                            constraint_map[key].lower_bounds.conservativeResize(constraint_map[key].lower_bounds.rows(), constraint_map[key].lower_bounds.cols() + constraint.lower_bounds.cols());
                            constraint_map[key].lower_bounds << constraint.lower_bounds;

                            constraint_map[key].upper_bounds.conservativeResize(constraint_map[key].upper_bounds.rows(), constraint_map[key].upper_bounds.cols() + constraint.upper_bounds.cols());
                            constraint_map[key].upper_bounds << constraint.upper_bounds;
                        }
                        else
                        {
                            // If the constraint doesn't exist, add it to the map
                            opt::constraint_evaluations_t new_constraint;
                            new_constraint.metadata.plot_groupings = {constraint.metadata.plot_groupings[i]};
                            new_constraint.metadata.plot_titles = {constraint.metadata.plot_titles[i]};
                            new_constraint.metadata.plot_names = {constraint.metadata.plot_names[i]};
                            new_constraint.metadata.name = constraint.metadata.name;

                            new_constraint.times = constraint.times;
                            new_constraint.evaluation = constraint.evaluation;
                            new_constraint.lower_bounds = constraint.lower_bounds;
                            new_constraint.upper_bounds = constraint.upper_bounds;

                            constraint_map[key] = new_constraint;
                        }
                    }
                }
            }

            std::cout << "Here" << std::endl;

            // Convert the map to a vector
            std::vector<opt::constraint_evaluations_t> new_constraints;
            for (auto &pair : constraint_map)
            {
                new_constraints.push_back(pair.second);
            }

            // Plot each group of constraints on the same graph
            for (size_t i = 0; i < new_constraints.size(); ++i)
            {
                std::vector<double> std_times(new_constraints[i].times.data(), new_constraints[i].times.data() + new_constraints[i].times.size());
                for (size_t j = 0; j < new_constraints[i].metadata.plot_titles.size(); ++j)
                {
                    gp << "set term wxt " << window_index << "\n"; // Create a new window for each plot
                    gp << "set xlabel 'Time'\n";
                    gp << "set ylabel 'Values'\n";
                    gp << "set title '" + new_constraints[i].metadata.plot_titles[j] + "'\n";
                    gp << "set style fill transparent solid 0.5 noborder\n";
                    std::string ss = "plot ";
                    for (size_t k = 0; k < new_constraints[i].metadata.plot_names[j].size(); ++k)
                    {
                        if (k != 0)
                        {
                            ss += ", ";
                        }
                        // Generate random values for the red, green, and blue components
                        int red = rand() % 256;
                        int green = rand() % 256;
                        int blue = rand() % 256;

                        // Convert the RGB values to hexadecimal
                        std::stringstream color;
                        color << "rgb '#" << std::setw(2) << std::setfill('0') << std::hex << red
                              << std::setw(2) << std::setfill('0') << std::hex << green
                              << std::setw(2) << std::setfill('0') << std::hex << blue << "'";

                        ss += "'-' with linespoints linestyle " + std::to_string(3 * k + 1) + " title '" + new_constraints[i].metadata.plot_names[j][k] + " Evaluation', ";
                        ss += "'-' with linespoints linestyle " + std::to_string(3 * k + 2) + " title '" + new_constraints[i].metadata.plot_names[j][k] + " Lower Bound', ";
                        ss += "'-' with linespoints linestyle " + std::to_string(3 * k + 3) + " title '" + new_constraints[i].metadata.plot_names[j][k] + " Upper Bound', ";
                        ss += "'-' with filledcurves lc " + color.str() + " title '" + new_constraints[i].metadata.plot_names[j][k] + " Bounds'";
                    }
                    ss += "\n";
                    gp << ss;
                    Eigen::MatrixXd block_eval = new_constraints[i].evaluation.block(0, std::get<0>(new_constraints[i].metadata.plot_groupings[j]), new_constraints[i].evaluation.rows(), std::get<1>(new_constraints[i].metadata.plot_groupings[j]) - std::get<0>(new_constraints[i].metadata.plot_groupings[j]));
                    Eigen::MatrixXd block_lb = new_constraints[i].lower_bounds.block(0, std::get<0>(new_constraints[i].metadata.plot_groupings[j]), new_constraints[i].lower_bounds.rows(), std::get<1>(new_constraints[i].metadata.plot_groupings[j]) - std::get<0>(new_constraints[i].metadata.plot_groupings[j]));
                    Eigen::MatrixXd block_ub = new_constraints[i].upper_bounds.block(0, std::get<0>(new_constraints[i].metadata.plot_groupings[j]), new_constraints[i].upper_bounds.rows(), std::get<1>(new_constraints[i].metadata.plot_groupings[j]) - std::get<0>(new_constraints[i].metadata.plot_groupings[j]));

                    for (Eigen::Index k = 0; k < block_eval.cols(); ++k)
                    {
                        Eigen::VectorXd colMatrix = block_eval.col(k).matrix();
                        std::vector<double> std_col_vector(colMatrix.data(), colMatrix.data() + colMatrix.size());
                        gp.send1d(std::make_tuple(std_times, std_col_vector));
                        colMatrix = block_lb.col(k).matrix();
                        std::vector<double> std_col_vector_lb(colMatrix.data(), colMatrix.data() + colMatrix.size());
                        gp.send1d(std::make_tuple(std_times, std_col_vector_lb));
                        colMatrix = block_ub.col(k).matrix();
                        std::vector<double> std_col_vector_ub(colMatrix.data(), colMatrix.data() + colMatrix.size());
                        gp.send1d(std::make_tuple(std_times, std_col_vector_ub));

                        // Send the data for the filled curves
                        gp.send1d(std::make_tuple(std_times, std_col_vector_lb, std_col_vector_ub));
                    }
                    ++window_index;
                }
            }
        }
    }
}