#pragma once

#include "third-party/gnuplot-iostream/gnuplot-iostream.h"
#include "galileo/opt/Solution.h"
#include <iostream>
#include <vector>

namespace galileo
{
    namespace tools
    {
        using tuple_size_t = std::tuple<size_t, size_t>;

        class GNUPlotInterface
        {
        public:
            GNUPlotInterface(opt::solution_t solution_, std::vector<std::vector<opt::constraint_evaluations_t>> constraint_evaluations_);

            void PlotSolution(std::vector<tuple_size_t> state_groups, std::vector<tuple_size_t> input_groups,
                              std::vector<std::string> state_title_names, std::vector<std::vector<std::string>> state_names,
                              std::vector<std::string> input_title_names, std::vector<std::vector<std::string>> input_names);

            void PlotConstraints();

        protected:
            opt::solution_t solution;
            std::vector<std::vector<opt::constraint_evaluations_t>> constraint_evaluations;

            Gnuplot gp;
            size_t window_index;
        };
    }
}