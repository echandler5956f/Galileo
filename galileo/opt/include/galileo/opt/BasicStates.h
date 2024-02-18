#pragma once

#include "galileo/opt/States.h"

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
             * @param nx Number of state variables
             * @param nu Number of control variables
             */
            BasicStates(int nx, int nu)
            {
                this->nx = nx;
                this->ndx = nx;
                this->nu = nu;
            }
        };
    }
}