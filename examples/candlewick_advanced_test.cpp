#include <candlewick/core/RenderContext.h>
#include <candlewick/core/GuiSystem.h>
#include <candlewick/core/DepthAndShadowPass.h>
#include <candlewick/core/LightUniforms.h>
#include <candlewick/core/CameraControls.h>
#include <candlewick/core/Components.h>

#include <candlewick/multibody/RobotScene.h>
#include <candlewick/multibody/RobotLoader.h>
#include <candlewick/primitives/Primitives.h>

#include <candlewick/utils/VideoRecorder.h>
#include <candlewick/utils/WriteTextureToImage.h>

#include <imgui.h>
#include <imgui_impl_sdl3.h>

#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/data.hpp>
#include <pinocchio/multibody/geometry.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/algorithm/geometry.hpp>

#include <SDL3/SDL_log.h>
#include <SDL3/SDL_init.h>
#include <SDL3/SDL_events.h>
#include <SDL3/SDL_assert.h>
#include <SDL3/SDL_gpu.h>

#include <CLI/App.hpp>
#include <CLI/Formatter.hpp>
#include <CLI/Config.hpp>

#include <entt/entity/registry.hpp>

namespace pin = pinocchio;
using namespace candlewick;
using multibody::RobotScene;
using multibody::RobotSpec;

const float kScrollZoom = 0.05f;

/// Application constants
constexpr Uint32 wWidth = 3840;
constexpr Uint32 wHeight = 2160;
constexpr float aspectRatio = float(wWidth) / float(wHeight);

/// Application state
static Radf currentFov = 55.0_degf;
static float nearZ = 0.01f;
static float farZ = 10.f;
static float currentOrthoScale = 1.f;
static CylindricalCamera g_camera{{
    .projection = perspectiveFromFov(currentFov, aspectRatio, nearZ, farZ),
    .view = Eigen::Isometry3f{lookAt({2.0, 0, 2.}, Float3::Zero())},
}};
static CameraProjection g_cameraType = CameraProjection::PERSPECTIVE;
static bool quitRequested = false;

static float pixelDensity;
static float displayScale;

static void updateFov(Radf newFov) {
  g_camera.camera.projection =
      perspectiveFromFov(newFov, aspectRatio, nearZ, farZ);
  currentFov = newFov;
}

static void updateOrtho(float zoom) {
  float iz = 1.f / zoom;
  g_camera.camera.projection =
      orthographicMatrix({iz * aspectRatio, iz}, -8., 8.);
  currentOrthoScale = zoom;
}

void eventLoop(const RenderContext &renderer) {
  pixelDensity = renderer.window.pixelDensity();
  displayScale = renderer.window.displayScale();
  const float rotSensitivity = 5e-3f * pixelDensity;
  const float panSensitivity = 1e-2f * pixelDensity;
  
  SDL_Event event;
  while (SDL_PollEvent(&event)) {
    ImGui_ImplSDL3_ProcessEvent(&event);
    ImGuiIO &io = ImGui::GetIO();
    
    if (event.type == SDL_EVENT_QUIT) {
      SDL_Log("Application exit requested.");
      quitRequested = true;
      break;
    }

    if (io.WantCaptureMouse | io.WantCaptureKeyboard)
      continue;
      
    switch (event.type) {
    case SDL_EVENT_MOUSE_WHEEL: {
      float wy = event.wheel.y;
      const float scaleFac = std::exp(kScrollZoom * wy);
      switch (g_cameraType) {
      case CameraProjection::ORTHOGRAPHIC:
        updateOrtho(std::clamp(scaleFac * currentOrthoScale, 0.1f, 2.f));
        break;
      case CameraProjection::PERSPECTIVE:
        updateFov(Radf(std::min(currentFov * scaleFac, Radf{170.0_degf})));
        break;
      }
      break;
    }
    case SDL_EVENT_KEY_DOWN: {
      const float step_size = 0.06f;
      switch (event.key.key) {
      case SDLK_LEFT:
        g_camera.localTranslate({+step_size, 0, 0});
        break;
      case SDLK_RIGHT:
        g_camera.localTranslate({-step_size, 0, 0});
        break;
      case SDLK_UP:
        g_camera.dolly(+step_size);
        break;
      case SDLK_DOWN:
        g_camera.dolly(-step_size);
        break;
      }
      break;
    }
    case SDL_EVENT_MOUSE_MOTION: {
      SDL_MouseButtonFlags mouseButton = event.motion.state;
      bool controlPressed = SDL_GetModState() & SDL_KMOD_CTRL;
      Float2 mvt{event.motion.xrel, event.motion.yrel};
      
      if (mouseButton & SDL_BUTTON_LMASK) {
        if (controlPressed) {
          g_camera.moveInOut(0.95f, event.motion.yrel);
        } else {
          g_camera.viewportDrag(mvt, rotSensitivity, panSensitivity);
        }
      }
      if (mouseButton & SDL_BUTTON_MMASK) {
        g_camera.pan(mvt, 5e-3f);
      }
      if (mouseButton & SDL_BUTTON_RMASK) {
        float camXLocRotSpeed = 0.01f * pixelDensity;
        camera_util::localRotateXAroundOrigin(g_camera, camXLocRotSpeed *
                                                            event.motion.yrel);
      }
      break;
    }
    }
  }
}

static void screenshot_button_callback(RenderContext &renderer,
                                       media::TransferBufferPool &pool,
                                       const char *filename) {
  const auto &device = renderer.device;
  CommandBuffer command_buffer{device};
  renderer.waitAndAcquireSwapchain(command_buffer);

  SDL_Log("Saving screenshot at %s", filename);
  media::saveTextureToFile(command_buffer, device, pool, renderer.swapchain,
                           renderer.getSwapchainTextureFormat(), wWidth,
                           wHeight, filename);
}

static const RobotSpec ur_robot_spec =
    RobotSpec{
        "urdf/ur5_gripper.urdf",
        "srdf/ur5_gripper.srdf",
        std::filesystem::path(EXAMPLE_ROBOT_DATA_MODEL_DIR).parent_path(),
        "robots/ur_description",
    }
        .ensure_absolute_filepaths();

int main(int argc, char **argv) {
  CLI::App app{"UR5 Advanced Rendering Example"};
  bool performRecording{false};
  
  argv = app.ensure_utf8(argv);
  app.add_flag("-r,--record", performRecording, "Record output");
  CLI11_PARSE(app, argc, argv);

  if (!SDL_Init(SDL_INIT_VIDEO))
    return 1;

  // Create render context with depth buffer
  RenderContext renderer{
      Device{auto_detect_shader_format_subset(), false},
      Window(__FILE__, wWidth, wHeight, 0),
      SDL_GPU_TEXTUREFORMAT_D16_UNORM,
  };

  entt::registry registry{};

  // Load robot model
  pin::Model model;
  pin::GeometryModel geom_model;
  loadModels(ur_robot_spec, model, &geom_model, NULL, true);

  pin::Data pin_data{model};
  pin::GeometryData geom_data{geom_model};

  // Configure advanced rendering features
  RobotScene::Config robot_scene_config;
  robot_scene_config.triangle_has_prepass = true;
  robot_scene_config.enable_normal_target = true;
  robot_scene_config.enable_msaa = true;
  robot_scene_config.msaa_samples = SDL_GPU_SAMPLECOUNT_8;

  // Create robot scene with advanced config
  RobotScene robot_scene{registry, renderer, geom_model, geom_data,
                         robot_scene_config};
  
  // Configure high-quality lighting
  robot_scene.directionalLight = {
      DirectionalLight{
          .direction = {-1.f, 0.f, -1.f},
          .color = {1.0, 1.0, 1.0},
          .intensity = 8.0,
      },
      DirectionalLight{
          .direction = {0.5, 1., -1.},
          .color = {1.0, 1.0, 1.0},
          .intensity = 8.0,
      },
  };

  // Add floor plane
  const Eigen::Affine3f plane_transform{Eigen::UniformScaling<float>(3.0f)};
  entt::entity plane_entity = robot_scene.addEnvironmentObject(
      loadPlaneTiled(0.5f, 20, 20), plane_transform.matrix());
  auto &plane_obj = registry.get<MeshMaterialComponent>(plane_entity);

  DepthPass depthPass(renderer.device, plane_obj.mesh.layout(),
                      renderer.depth_texture,
                      {SDL_GPU_CULLMODE_NONE, 0.05f, 0.f, true, false});

  const char *screenshot_filename = nullptr;

  // GUI system
  GuiSystem gui_system{
      renderer, [&](const RenderContext &r) {
        IMGUI_CHECKVERSION();

        static bool show_about_window = false;
        static bool show_imgui_window = false;

        if (show_about_window)
          showCandlewickAboutWindow(&show_about_window);
        if (show_imgui_window)
          ImGui::ShowAboutWindow(&show_imgui_window);

        ImGuiWindowFlags window_flags = ImGuiWindowFlags_AlwaysAutoResize | 
                                       ImGuiWindowFlags_MenuBar;
        ImGui::SetNextWindowPos({20, 20}, ImGuiCond_FirstUseEver);
        ImGui::Begin("Advanced Rendering Controls", nullptr, window_flags);

        if (ImGui::BeginMenuBar()) {
          ImGui::MenuItem("About Dear ImGui", NULL, &show_imgui_window);
          ImGui::MenuItem("About Candlewick", NULL, &show_about_window);
          ImGui::EndMenuBar();
        }

        ImGui::Text("Video driver: %s", SDL_GetCurrentVideoDriver());
        ImGui::SameLine();
        ImGui::Text("Device driver: %s", r.device.driverName());
        ImGui::Text("Display pixel density: %.2f / scale: %.2f",
                    r.window.pixelDensity(), r.window.displayScale());

        ImGui::SeparatorText("Camera");
        bool ortho_change, persp_change;
        ortho_change = ImGui::RadioButton("Orthographic", (int *)&g_cameraType,
                                          int(CameraProjection::ORTHOGRAPHIC));
        ImGui::SameLine();
        persp_change = ImGui::RadioButton("Perspective", (int *)&g_cameraType,
                                          int(CameraProjection::PERSPECTIVE));
        
        switch (g_cameraType) {
        case CameraProjection::ORTHOGRAPHIC:
          ortho_change |=
              ImGui::DragFloat("zoom", &currentOrthoScale, 0.01f, 0.1f, 2.f,
                               "%.3f", ImGuiSliderFlags_AlwaysClamp);
          if (ortho_change)
            updateOrtho(currentOrthoScale);
          break;
        case CameraProjection::PERSPECTIVE:
          Degf newFov{currentFov};
          persp_change |=
              ImGui::DragFloat("fov", newFov, 1.f, 15.f, 90.f, "%.3f",
                               ImGuiSliderFlags_AlwaysClamp);
          persp_change |=
              ImGui::SliderFloat("Near plane", &nearZ, 0.01f, 0.8f * farZ);
          persp_change |= ImGui::SliderFloat("Far plane", &farZ, nearZ, 20.f);
          if (persp_change)
            updateFov(Radf(newFov));
          break;
        }

        ImGui::SeparatorText("Rendering Features");
        ImGui::Checkbox("Ambient occlusion (SSAO)", &robot_scene.config().enable_ssao);
        ImGui::Checkbox("Depth pre-pass", &robot_scene.config().triangle_has_prepass);

        ImGui::SeparatorText("Lights");
        guiAddLightControls(robot_scene.directionalLight, robot_scene.numLights());

        ImGui::SeparatorText("Materials");
        ImGui::ColorEdit4("Floor color", plane_obj.materials[0].baseColor.data());

        ImGui::SeparatorText("Screenshots");
        static std::string scr_filename;
        guiAddFileDialog(renderer.window, DialogFileType::IMAGES, scr_filename);
        if (ImGui::Button("Take screenshot")) {
          if (scr_filename.empty())
            generateMediaFilenameFromTimestamp("ur5_advanced", scr_filename);
          screenshot_filename = scr_filename.c_str();
        }

        ImGui::End();
      }};

  // Main application loop
  Uint32 frameNo = 0;
  std::srand(42);
  Eigen::VectorXd q0 = pin::neutral(model);
  Eigen::VectorXd q1 = pin::randomConfiguration(model);

  media::TransferBufferPool transfer_buffer_pool{renderer.device};
  media::VideoRecorder recorder{NoInit};
  if (performRecording) {
    media::VideoRecorder::Settings settings;
    settings.fps = 50;
    recorder.open(wWidth, wHeight, "ur5_advanced.mp4", settings);
  }

  Eigen::VectorXd q = q0;
  Eigen::VectorXd qn = q;
  Eigen::VectorXd v{model.nv};
  const double dt = 1e-2;

  while (!quitRequested) {
    // Animation logic
    eventLoop(renderer);
    double alpha = 0.5 * (1. + std::sin(frameNo * dt));
    pin::interpolate(model, q0, q1, alpha, qn);
    v = pin::difference(model, q, qn) / dt;
    pin::forwardKinematics(model, pin_data, qn, v);
    pin::updateFramePlacements(model, pin_data);
    pin::updateGeometryPlacements(model, pin_data, geom_model, geom_data);
    q = qn;
    robot_scene.updateTransforms();

    // Rendering
    CommandBuffer command_buffer = renderer.acquireCommandBuffer();

    if (renderer.waitAndAcquireSwapchain(command_buffer)) {
      const GpuMat4 viewProj = g_camera.camera.viewProj();
      
      // Collect shadow casters
      robot_scene.collectOpaqueCastables();
      auto &castables = robot_scene.castables();
      
      renderShadowPassFromFrustum(command_buffer, robot_scene.shadowPass,
                                  robot_scene.directionalLight, castables,
                                  frustumFromCameraViewProj(viewProj));
      
      depthPass.render(command_buffer, viewProj, castables);
      
      robot_scene.renderOpaque(command_buffer, g_camera);
      robot_scene.renderTransparent(command_buffer, g_camera);
      
      gui_system.render(command_buffer);
    } else {
      SDL_Log("Failed to acquire swapchain: %s", SDL_GetError());
      continue;
    }

    command_buffer.submit();

    // Video recording
    if (performRecording) {
      CommandBuffer command_buffer = renderer.acquireCommandBuffer();
      auto swapchain_format = renderer.getSwapchainTextureFormat();
      recorder.writeTextureToVideoFrame(command_buffer, renderer.device,
                                        transfer_buffer_pool,
                                        renderer.swapchain, swapchain_format);
    }
    
    // Screenshots
    if (screenshot_filename) {
      screenshot_button_callback(renderer, transfer_buffer_pool,
                                 screenshot_filename);
      screenshot_filename = nullptr;
    }
    
    frameNo++;
  }

  // Cleanup
  SDL_WaitForGPUIdle(renderer.device);
  depthPass.release();  // Clean up depth pass
  robot_scene.release();
  gui_system.release();
  transfer_buffer_pool.release();
  recorder.close();
  renderer.destroy();
  SDL_Quit();
  return 0;
} 