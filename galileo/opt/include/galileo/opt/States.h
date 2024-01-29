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

            using FactoryFunction = std::function<std::unique_ptr<States>(std::vector<int>)>;

            static void registerRobotType(const std::string &type, FactoryFunction function)
            {
                getFactories()[type] = function;
            }
            static std::unique_ptr<States> create(const std::string &type, std::vector<int> args)
            {
                return getFactories()[type](args);
            }

        private:
            static std::map<std::string, FactoryFunction> &getFactories()
            {
                static std::map<std::string, FactoryFunction> factories;
                return factories;
            }
        };
    }
}