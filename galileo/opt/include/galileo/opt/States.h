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

            static void registerRobotType(const std::string &type, FactoryFunction function);
            static std::unique_ptr<States> create(const std::string &type, std::vector<int> args);

        private:
            static std::map<std::string, FactoryFunction> factories;
        };

        std::map<std::string, States::FactoryFunction> States::factories;

        void States::registerRobotType(const std::string &type, FactoryFunction function)
        {
            factories[type] = function;
        }

        std::unique_ptr<States> States::create(const std::string &type, std::vector<int> args)
        {
            return factories[type](args);
        }
    }
}