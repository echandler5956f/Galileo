#pragma once

#include "third-party/gnuplot-iostream/gnuplot-iostream.h"
#include "galileo/opt/Solution.h"
#include <iostream>
#include <vector>
#include <cstdlib> // for rand()
#include <sstream> // for std::stringstream
#include <iomanip> // for std::setw and std::setfill

namespace galileo
{
    namespace tools
    {
        using tuple_size_t = std::tuple<size_t, size_t>;

        /**
         * @brief Class to interface with GNUPlot, used to plot solutions and constraints.
         * 
         */
        class GNUPlotInterface
        {
        public:
            /**
             * @brief Construct a new GNUPlotInterface object.
             * 
             * @param dir_ Directory where the plots will be saved.
             */
            GNUPlotInterface(std::string dir_);

            /**
             * @brief Plot the solution.
             * 
             * @param solution The solution to be plotted.
             * @param state_groups A set of ranges indicating the groups of states to be plotted.
             * @param input_groups A set of ranges indicating the groups of inputs to be plotted.
             * @param state_title_names The titles of the states groups. Used for the title.
             * @param state_names The names of the states. Used for the legend.
             * @param input_title_names The titles of the inputs groups. Used for the title.
             * @param input_names The names of the inputs. Used for the legend.
             */
            void PlotSolution(opt::solution::solution_t solution,
                              std::vector<tuple_size_t> state_groups, std::vector<tuple_size_t> input_groups,
                              std::vector<std::string> state_title_names, std::vector<std::vector<std::string>> state_names,
                              std::vector<std::string> input_title_names, std::vector<std::vector<std::string>> input_names);

            /**
             * @brief Plot the constraints.
             * 
             * @param constraint_evaluations The constraint evaluations to be plotted.
             */
            void PlotConstraints(std::vector<std::vector<opt::constraint_evaluations_t>> constraint_evaluations);

        protected:
            /**
             * @brief The Gnuplot object.
             * 
             */
            Gnuplot gp_;

            /**
             * @brief The index of the current plotting window.
             * 
             */
            size_t window_index_ = 0;

            /**
             * @brief The directory where the plots will be saved.
             * 
             */
            std::string dir_;
        };
    }
}