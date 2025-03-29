#ifndef __galileo_multibody_robot_robot_data_base_hpp__
#define __galileo_multibody_robot_robot_data_base_hpp__

#include "galileo/multibody/robot/robot-base.hpp"
#include "galileo/multibody/robot/robot-model-base.hpp"

namespace galileo
{
    namespace multibody
    {

        template <typename Derived>
        struct RobotDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using RobotDerived = typename traits<Derived>::RobotDerived;
            GALILEO_ROBOT_BASIC_TYPEDEF(RobotDerived);
            GALILEO_ROBOT_CONSTANTS(RobotDerived);
            GALILEO_ROBOT_DATA_TYPEDEF(RobotDerived);

        protected:
            inline RobotDataBase()
            {
            }

            inline RobotDataBase(const RobotDataBase &clone)
            {
                *this = clone;
            }

            inline RobotDataBase &operator=(const RobotDataBase &clone)
            {
                return *this;
            }

        }; // struct RobotDataBase

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_robot_robot_data_base_hpp__
