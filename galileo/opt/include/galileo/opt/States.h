#pragma once

#include <memory>
#include <map>
#include <functional>
#include <vector>

using namespace std;

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

            using FactoryFunction = function<unique_ptr<States>(vector<int>)>;

            static void registerRobotType(const string &type, FactoryFunction function)
            {
                getFactories()[type] = function;
            }
            static unique_ptr<States> create(const string &type, vector<int> args)
            {
                return getFactories()[type](args);
            }

        private:
            static map<string, FactoryFunction> &getFactories()
            {
                static map<string, FactoryFunction> factories;
                return factories;
            }
        };
    }
}