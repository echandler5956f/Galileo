#include "candlewick/multibody/Visualizer.h"
#include "candlewick/multibody/RobotLoader.h"

#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/algorithm/geometry.hpp>
#include <chrono>

#include <CLI/App.hpp>
#include <CLI/Formatter.hpp>
#include <CLI/Config.hpp>

using namespace candlewick::multibody;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

static const RobotSpec go1_robot_spec =
    RobotSpec{
        "urdf/go1.urdf",
        "srdf/go1.srdf",
        fs::path(EXAMPLE_ROBOT_DATA_MODEL_DIR).parent_path(),
        "robots/go1_description",
        true
    }
        .ensure_absolute_filepaths();

int main(int argc, char **argv) {
  CLI::App app{"Visualizer example"};
  argv = app.ensure_utf8(argv);
  std::array<Uint32, 2> window_dims{1920u, 1080u};
  double fps;

  app.add_option("--dims", window_dims, "Window dimensions.")
      ->capture_default_str();
  app.add_option<double, unsigned int>("--fps", fps, "Framerate")
      ->default_val(60);

  CLI11_PARSE(app, argc, argv);

  pinocchio::Model model;
  pinocchio::GeometryModel geom_model;
  loadModels(go1_robot_spec, model, &geom_model, NULL);

  Visualizer visualizer{{window_dims[0], window_dims[1]}, model, geom_model};
  assert(!visualizer.hasExternalData());
  pinocchio::Data &data = visualizer.data();

  Eigen::VectorXd q0 = model.referenceConfigurations["standing"];
  Eigen::VectorXd q1 = pinocchio::neutral(model);

  std::cout << "q0: " << q0.transpose() << std::endl;
  std::cout << "q1: " << q1.transpose() << std::endl;

  double dt = 1. / static_cast<double>(fps);
  using duration_t = std::chrono::duration<double>;
  Eigen::VectorXd q = q0;
  Eigen::VectorXd qn = q;
  Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);

  double t = 0.;

  while (!visualizer.shouldExit()) {
    const auto now = steady_clock::now();

    double alpha = 0.5 * (std::sin(t) + 1.0);
    pinocchio::interpolate(model, q0, q1, alpha, q);

    pinocchio::forwardKinematics(model, data, q);
    auto base_to_rl_foot = data.oMf[model.getFrameId("base")].inverse() * data.oMf[model.getFrameId("RL_foot")];
    q[2] = -base_to_rl_foot.translation()[2];

    pinocchio::difference(model, qn, q, v);
    v /= dt;
    pinocchio::forwardKinematics(model, data, q, v);

    visualizer.display();
    std::this_thread::sleep_until(now + duration_t(dt));

    t += dt;
    qn = q;
  }
  return 0;
}
