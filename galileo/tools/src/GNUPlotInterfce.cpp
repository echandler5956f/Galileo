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

            Gnuplot gp;  // Create a single Gnuplot process here

            for (size_t i = 0; i < state_groups.size(); ++i)
            {
                gp << "set term wxt " << i << "\n";  // Create a new window for each plot
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
            }

            for (size_t i = 0; i < input_groups.size(); ++i)
            {
                gp << "set term wxt " << (i + state_groups.size()) << "\n";  // Create a new window for each plot
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
            }
        }

        void GNUPlotInterface::PlotConstraints()
        {
            std::vector<opt::constraint_evaluations_t> new_constraints;
            std::map<std::string, opt::constraint_evaluations_t *> constraint_map;

            for (auto &cons_in_phase : constraint_evaluations)
            {
                for (auto &constraint : cons_in_phase)
                {
                    if (constraint_map.count(constraint.name) == 0)
                    {
                        // This is a new constraint, add it to new_constraints and constraint_map
                        new_constraints.push_back(constraint);
                        constraint_map[constraint.name] = &new_constraints.back();
                    }
                    else
                    {
                        // This is an existing constraint, append the times and evaluation_and_bounds data
                        opt::constraint_evaluations_t *existing_constraint = constraint_map[constraint.name];
                        existing_constraint->times.conservativeResize(existing_constraint->times.size() + constraint.times.size());
                        existing_constraint->times.tail(constraint.times.size()) = constraint.times;
                        existing_constraint->evaluation_and_bounds.conservativeResize(existing_constraint->evaluation_and_bounds.rows(), existing_constraint->evaluation_and_bounds.cols() + constraint.evaluation_and_bounds.cols());
                        existing_constraint->evaluation_and_bounds.rightCols(constraint.evaluation_and_bounds.cols()) = constraint.evaluation_and_bounds;
                    }
                }
            }
            std::cout << "New constraints size: " << new_constraints.size() << std::endl;
            // Plot each group of constraints on the same graph
            for (size_t i = 0; i < new_constraints.size(); ++i)
            {
                std::cout << "Plotting constraint: " << new_constraints[i].name << std::endl;
                std::vector<double> std_times(new_constraints[i].times.data(), new_constraints[i].times.data() + new_constraints[i].times.size());
                Gnuplot gp;
                gp << "set xlabel 'Time'\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + new_constraints[i].name + "'\n";
                std::string ss = "plot ";
                std::string tmp_name;
                for (size_t j = 0; j < 3; ++j)
                {
                    if (j == 0)
                        tmp_name = "Evaluation";
                    else if (j == 1)
                        tmp_name = "Lower Bound";
                    else if (j == 2)
                        tmp_name = "Upper Bound";

                    if (j != 0)
                    {
                        ss += ", ";
                    }
                    ss += "'-' with linespoints linestyle " + std::to_string(j + 1) + " title '" + tmp_name + "'";
                }
                ss += "\n";
                gp << ss;
                for (Eigen::Index j = 0; j < new_constraints[i].evaluation_and_bounds.rows(); ++j)
                {
                    Eigen::MatrixXd rowMatrix = new_constraints[i].evaluation_and_bounds.row(j).matrix();
                    std::vector<double> std_row_vector(rowMatrix.data(), rowMatrix.data() + rowMatrix.size());
                    gp.send1d(std::make_tuple(std_times, std_row_vector));
                }
            }
        }
    }
}