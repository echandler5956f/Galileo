#pragma once

#include <memory>
#include <map>
#include <vector>

namespace galileo
{
    namespace opt
    {
        /**
         * @brief Simple slicer class for getting state variables.
         *
         */
        class States
        {
        public:
            /**
             * @brief Destroy the States object.
             * 
             */
            virtual ~States() = default;

            /**
             * @brief Number of state variables.
             *
             */
            int nx;

            /**
             * @brief Number of state deviant variables.
             *
             */
            int ndx;

            /**
             * @brief Number of input variables.
             *
             */
            int nu;
        };

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
            BasicStates(int nx_, int nu_)
            {
                this->nx = nx_;
                this->ndx = nx_;
                this->nu = nu_;
            }
        };
    }
}