#include "galileo/tools/GNUPlotInterface.h"

namespace galileo
{
    namespace tools
    {
        GNUPlotInterface::GNUPlotInterface(opt::solution_t solution_, std::vector<std::vector<opt::constraint_evaluations_t>> constraint_evaluations_)
        {
            this->solution = solution_;
            this->constraint_evaluations = constraint_evaluations_;
            std::cout << "Solution state size: " << solution.state_result.rows() << "x" << solution.state_result.cols() << std::endl;
            std::cout << "Solution input size: " << solution.input_result.rows() << "x" << solution.input_result.cols() << std::endl;
            std::cout << "Solution time size: " << solution.times.size() << std::endl;
        }

        void GNUPlotInterface::PlotSolution(std::vector<tuple_size_t> state_groups, std::vector<tuple_size_t> input_groups,
                                            std::vector<std::string> state_title_names, std::vector<std::vector<std::string>> state_names,
                                            std::vector<std::string> input_title_names, std::vector<std::vector<std::string>> input_names)
        {
            std::vector<double> std_times(solution.times.data(), solution.times.data() + solution.times.size());
            for (size_t i = 0; i < state_groups.size(); ++i)
            {
                Gnuplot gp;
                gp << "set xlabel 'Time'\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + state_title_names[i] + "'\n";
                std::string ss = "plot";
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
                Eigen::MatrixXd block = solution.state_result.block(std::get<0>(state_groups[i]), 0, std::get<1>(state_groups[i]), solution.state_result.cols());
                for (Eigen::Index j = 0; j < block.rows(); ++j)
                {
                    Eigen::MatrixXd rowMatrix = block.row(j).matrix();
                    std::vector<double> std_row_vector(rowMatrix.data(), rowMatrix.data() + rowMatrix.size());
                    gp.send1d(std::make_tuple(std_times, std_row_vector));
                }
            }
            for (size_t i = 0; i < input_groups.size(); ++i)
            {
                Gnuplot gp;
                gp << "set xlabel 'Time'\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + input_title_names[i] + "'\n";
                std::string ss = "plot";
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
                Eigen::MatrixXd block = solution.input_result.block(std::get<0>(input_groups[i]), 0, std::get<1>(input_groups[i]), solution.input_result.cols());
                for (Eigen::Index j = 0; j < block.rows(); ++j)
                {
                    Eigen::MatrixXd rowMatrix = block.row(j).matrix();
                    std::vector<double> std_row_vector(rowMatrix.data(), rowMatrix.data() + rowMatrix.size());
                    gp.send1d(std::make_tuple(std_times, std_row_vector));
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

            // Plot each group of constraints on the same graph
            for (size_t i = 0; i < new_constraints.size(); ++i)
            {
                std::vector<double> std_times(new_constraints[i].times.data(), new_constraints[i].times.data() + new_constraints[i].times.size());
                Gnuplot gp;
                gp << "set xlabel 'Time'\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + new_constraints[i].name + "'\n";
                std::string ss = "plot";
                std::string tmp_name ;
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