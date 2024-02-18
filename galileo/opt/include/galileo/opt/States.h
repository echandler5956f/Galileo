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
    }
}