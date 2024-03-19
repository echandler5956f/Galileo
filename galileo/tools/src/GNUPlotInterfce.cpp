#include "galileo/tools/GNUPlotInterface.h"

namespace galileo
{
    namespace tools
    {
        GNUPlotInterface::GNUPlotInterface(opt::solution::solution_t solution_, std::vector<std::vector<opt::solution::constraint_evaluations_t>> constraint_evaluations_)
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
            // gp << "set terminal qt " << 0 << "\n"; // Create a dummy window
            // gp << "plot '-' with lines linestyle 1 title 'Dummy'\n";
            // gp.send1d(std::make_tuple(std_times, std_times));

            window_index = 0;
            for (size_t i = 0; i < state_groups.size(); ++i)
            {
                std::string filename = "output_" + state_title_names[i] + ".png";
                gp << "set terminal pngcairo\n";
                gp << "set output '" << filename << "'\n";
                // gp << "set terminal qt " << window_index << "\n"; // Create a new window for each plot
                gp << "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#FFFFFF' behind\n";
                gp << "set xlabel 'Time'\n";
                gp << "set xrange [" << solution.times(0) << ":" << solution.times(solution.times.size() - 1) << "]\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + state_title_names[i] + "'\n";
                // gp << "set key outside\n"; // Move the legend outside of the plot
                std::string ss = "plot ";
                for (size_t j = 0; j < state_names[i].size(); ++j)
                {
                    if (j != 0)
                    {
                        ss += ", ";
                    }
                    ss += "'-' with lines linestyle " + std::to_string(j + 1) + " title '" + state_names[i][j] + "'";
                }
                ss += "\n";
                gp << ss;
                sleep(0.1);
                Eigen::MatrixXd block = solution.state_result.block(0, std::get<0>(state_groups[i]), solution.state_result.rows(), std::get<1>(state_groups[i]) - std::get<0>(state_groups[i]));
                for (Eigen::Index j = 0; j < block.cols(); ++j)
                {
                    Eigen::VectorXd colMatrix = block.col(j).matrix();
                    std::vector<double> std_col_vector(colMatrix.data(), colMatrix.data() + colMatrix.size());
                    gp.send1d(std::make_tuple(std_times, std_col_vector));
                }
                // gp << "set output\n"; // Close the output file
                ++window_index;
            }

            for (size_t i = 0; i < input_groups.size(); ++i)
            {
                std::string filename = "output_" + input_title_names[i] + ".png";
                gp << "set terminal pngcairo\n";
                gp << "set output '" << filename << "'\n";
                // gp << "set terminal qt " << window_index << "\n"; // Create a new window for each plot
                gp << "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#FFFFFF' behind\n";
                gp << "set xlabel 'Time'\n";
                gp << "set xrange [" << solution.times(0) << ":" << solution.times(solution.times.size() - 1) << "]\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + input_title_names[i] + "'\n";
                // gp << "set key outside\n"; // Move the legend outside of the plot
                std::string ss = "plot ";
                for (size_t j = 0; j < input_names[i].size(); ++j)
                {
                    if (j != 0)
                    {
                        ss += ", ";
                    }
                    ss += "'-' with lines linestyle " + std::to_string(j + 1) + " title '" + input_names[i][j] + "'";
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
                // gp << "set output\n"; // Close the output file
                ++window_index;
            }
        }

        void GNUPlotInterface::PlotConstraints()
        {
            std::map<std::string, std::vector<opt::solution::constraint_evaluations_t>> constraint_map;

            // For each phase
            for (auto &cons_in_phase : constraint_evaluations)
            {
                // For each constraint type in the phase
                for (auto &constraint_type : cons_in_phase)
                {
                    // Constraints in different phases can have the same name; we want to combine them according to their plot_groupings
                    // For each sub constraint in the constraint type
                    for (size_t i = 0; i < constraint_type.metadata.plot_titles.size(); ++i)
                    {
                        std::string key = constraint_type.metadata.name + constraint_type.metadata.plot_titles[i];
                        size_t s1 = std::get<0>(constraint_type.metadata.plot_groupings[i]);
                        size_t s2 = std::get<1>(constraint_type.metadata.plot_groupings[i]);

                        opt::solution::constraint_evaluations_t new_constraint;
                        new_constraint.metadata.plot_groupings = {constraint_type.metadata.plot_groupings[i]};
                        new_constraint.metadata.plot_titles = {constraint_type.metadata.plot_titles[i]};
                        new_constraint.metadata.plot_names = {constraint_type.metadata.plot_names[i]};
                        new_constraint.metadata.name = constraint_type.metadata.name;
                        new_constraint.times = constraint_type.times;
                        new_constraint.evaluation = constraint_type.evaluation.block(0, s1, constraint_type.evaluation.rows(), s2 - s1);
                        new_constraint.lower_bounds = constraint_type.lower_bounds.block(0, s1, constraint_type.lower_bounds.rows(), s2 - s1);
                        new_constraint.upper_bounds = constraint_type.upper_bounds.block(0, s1, constraint_type.upper_bounds.rows(), s2 - s1);
                        if (constraint_map.count(key) > 0)
                        {
                            constraint_map[key].push_back(new_constraint);
                        }
                        else
                        {
                            constraint_map[key] = {new_constraint};
                        }
                    }
                }
            }

            // Convert the map to a vector
            std::vector<std::vector<opt::solution::constraint_evaluations_t>> new_constraints;
            for (auto &pair : constraint_map)
            {
                new_constraints.push_back(pair.second);
            }

            std::vector<std::string> colors;
            // Plot each group of constraints on the same graph
            for (size_t i = 0; i < new_constraints.size(); ++i)
            {
                std::string filename = "output_" + new_constraints[i][0].metadata.plot_titles[0] + ".png";
                gp << "set terminal pngcairo\n";
                gp << "set output '" << filename << "'\n";
                // gp << "set terminal qt " << window_index << "\n"; // Create a new window for each plot
                gp << "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#FFFFFF' behind\n";
                gp << "set xlabel 'Time'\n";
                gp << "set xrange [" << solution.times(0) << ":" << solution.times(solution.times.size() - 1) << "]\n";
                gp << "set ylabel 'Values'\n";
                gp << "set title '" + new_constraints[i][0].metadata.plot_titles[0] + "'\n";
                gp << "set style fill transparent solid 0.5 noborder\n";
                gp << "set key inside left\n"; // Move the legend outside of the plot
                std::string ss = "plot ";

                colors.clear();
                size_t plt_size = new_constraints[i][0].metadata.plot_names[0].size();
                size_t num_colors = 4 * plt_size;
                for (size_t j = 0; j < num_colors; ++j)
                {
                    // Generate random values for the red, green, and blue components
                    int red = rand() % 256;
                    int green = rand() % 256;
                    int blue = rand() % 256;

                    // Convert the RGB values to hexadecimal
                    std::stringstream color;
                    color << " lc rgb '#" << std::setw(2) << std::setfill('0') << std::hex << red
                          << std::setw(2) << std::setfill('0') << std::hex << green
                          << std::setw(2) << std::setfill('0') << std::hex << blue << "'";
                    colors.push_back(color.str());
                }

                for (size_t j = 0; j < new_constraints[i].size(); ++j)
                {
                    for (size_t k = 0; k < plt_size; ++k)
                    {
                        size_t idx = j * plt_size;
                        if (idx + k != 0)
                        {
                            ss += ", ";
                        }

                        if (j == 0)
                        {
                            ss += "'-' with lines linestyle " + std::to_string(idx + 3 * k + 1) + colors[4 * k + 0] + " title '" + new_constraints[i][0].metadata.plot_names[0][k] + " Evaluation', ";
                            ss += "'-' with lines linestyle " + std::to_string(idx + 3 * k + 2) + colors[4 * k + 1] + " title '" + new_constraints[i][0].metadata.plot_names[0][k] + " Lower Bound', ";
                            ss += "'-' with lines linestyle " + std::to_string(idx + 3 * k + 3) + colors[4 * k + 2] + " title '" + new_constraints[i][0].metadata.plot_names[0][k] + " Upper Bound', ";
                            ss += "'-' with filledcurves" + colors[4 * k + 3] + " title '" + new_constraints[i][0].metadata.plot_names[0][k] + " Bounds'";
                        }
                        else
                        {
                            ss += "'-' with lines linestyle " + std::to_string(idx + 3 * k + 1) + colors[4 * k + 0] + " notitle, ";
                            ss += "'-' with lines linestyle " + std::to_string(idx + 3 * k + 2) + colors[4 * k + 1] + " notitle, ";
                            ss += "'-' with lines linestyle " + std::to_string(idx + 3 * k + 3) + colors[4 * k + 2] + " notitle, ";
                            ss += "'-' with filledcurves" + colors[4 * k + 3] + " notitle";
                        }
                    }
                }
                ss += "\n";
                gp << ss;
                for (size_t j = 0; j < new_constraints[i].size(); ++j)
                {
                    std::vector<double> std_times(new_constraints[i][j].times.data(), new_constraints[i][j].times.data() + new_constraints[i][j].times.size());

                    Eigen::MatrixXd eval = new_constraints[i][j].evaluation;
                    Eigen::MatrixXd lb = new_constraints[i][j].lower_bounds;
                    Eigen::MatrixXd ub = new_constraints[i][j].upper_bounds;

                    for (Eigen::Index k = 0; k < eval.cols(); ++k)
                    {
                        Eigen::VectorXd colMatrix = eval.col(k).matrix();
                        std::vector<double> std_col_vector(colMatrix.data(), colMatrix.data() + colMatrix.size());
                        gp.send1d(std::make_tuple(std_times, std_col_vector));
                        colMatrix = lb.col(k).matrix();
                        std::vector<double> std_col_vector_lb(colMatrix.data(), colMatrix.data() + colMatrix.size());
                        gp.send1d(std::make_tuple(std_times, std_col_vector_lb));
                        colMatrix = ub.col(k).matrix();
                        std::vector<double> std_col_vector_ub(colMatrix.data(), colMatrix.data() + colMatrix.size());
                        gp.send1d(std::make_tuple(std_times, std_col_vector_ub));
                        // Send the data for the filled curves
                        gp.send1d(std::make_tuple(std_times, std_col_vector_lb, std_col_vector_ub));
                        double constraint_violation = 0.;
                        for (Eigen::Index cnt = 0; cnt < eval.rows(); ++cnt)
                        {
                            constraint_violation += pow(std::max(0., std_col_vector[cnt] - std_col_vector_ub[cnt]), 2);
                            constraint_violation += pow(std::max(0., std_col_vector_lb[cnt] - std_col_vector[cnt]), 2);
                        }
                        constraint_violation = sqrt(constraint_violation / eval.rows());
                        std::cout << "Constraint violation of " << new_constraints[i][0].metadata.plot_titles[0] << " " << new_constraints[i][j].metadata.name << ": " << constraint_violation << std::endl;
                    }
                }
                ++window_index;
            }
        }
    }
}