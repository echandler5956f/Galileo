#pragma once
#include "galileo/opt/States.h"

using namespace casadi;
using namespace std;

namespace galileo
{
    namespace opt
    {
        /**
         * @brief Simple slicer class for getting state variables.
         *
         */
        class BasicStates : public States
        {
        public:
            /**
             * @brief Construct a new Test State object
             *
             * @param args The number of position and velocity variables.
             */
            BasicStates(vector<int> args)
            {
                this->nx = args[0];
                this->ndx = args[0];
                this->nu = args[1];
            }
        };
        // Register the factory function for LeggedRobotStates
        bool registered = []
        {
            States::registerRobotType("BasicStates", [](std::vector<int> args)
                                      { return std::make_unique<BasicStates>(args); });
            return true;
        }();
    }
}