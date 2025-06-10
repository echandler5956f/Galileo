#ifndef __galileo_tests_models_pendulum_hpp__
#define __galileo_tests_models_pendulum_hpp__

#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>
#include <string>

namespace galileo {
namespace tests {

inline void build_pendulum(pinocchio::Model& model) {
    const std::string urdf_model = R"(
<robot name="simple_pendulum">
  <link name="base">
    <inertial>
      <mass value="1.0"/>
      <inertia ixx="1.0" ixy="0.0" ixz="0.0" iyy="1.0" iyz="0.0" izz="1.0"/>
    </inertial>
  </link>
  <link name="link1">
    <inertial>
      <mass value="1.0"/>
      <origin xyz="0.0 0.0 -0.5"/>
      <inertia ixx="0.25" ixy="0.0" ixz="0.0" iyy="0.25" iyz="0.0" izz="0.001"/>
    </inertial>
  </link>
  <joint name="joint1" type="revolute">
    <parent link="base"/>
    <child link="link1"/>
    <axis xyz="1 0 0"/>
    <limit effort="10" velocity="10"/>
  </joint>
</robot>
)";
    pinocchio::urdf::buildModelFromXML(urdf_model, model);
}

} // namespace tests
} // namespace galileo

#endif // __galileo_tests_models_pendulum_hpp__ 