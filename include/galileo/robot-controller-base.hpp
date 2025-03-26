#ifndef __galileo_robot_controller_base_hpp__
#define __galileo_robot_controller_base_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    template <typename Derived>
    class RobotControllerBase : internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RobotControllerDerived = typename traits<Derived>::RobotControllerDerived;
        GALILEO_ROBOT_CONTROLLER_BASIC_TYPEDEF(RobotControllerDerived);

        bool init();
        void update(const NumScalar &time);
        void starting(const NumScalar &time);
        void stopping(const NumScalar &time);

        void setupStateEstimation();
        void updateStateEstimation(const NumScalar &time);

        void setupMPC();
        void updateMPC(const NumScalar &time);

        void setupWBC();
        void updateWBC(const NumScalar &time);

        void setupSimulator();
        void updateSimulator(const NumScalar &time);

        void setupVisualizer();
        void updateVisualizer(const NumScalar &time);

    }; // class RobotControllerBase

} // namespace galileo

#endif // __galileo_robot_controller_base_hpp__