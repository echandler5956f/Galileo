#pragma once

#include <memory>
#include <map>
#include <functional>
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