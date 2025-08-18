set(IMGUI_DIR ${PROJECT_SOURCE_DIR}/external/imgui)
set(
  IMGUI_SRC
  ${IMGUI_DIR}/imgui.cpp
  ${IMGUI_DIR}/imgui_demo.cpp
  ${IMGUI_DIR}/imgui_draw.cpp
  ${IMGUI_DIR}/imgui_tables.cpp
  ${IMGUI_DIR}/imgui_widgets.cpp
  ${IMGUI_DIR}/backends/imgui_impl_sdl3.cpp
  ${IMGUI_DIR}/backends/imgui_impl_sdlgpu3.cpp
  PARENT_SCOPE
)

set(IMGUI_HEADERS ${IMGUI_DIR}/imconfig.h ${IMGUI_DIR}/imgui.h)

add_library(imgui_headers INTERFACE)
target_link_libraries(imgui_headers INTERFACE SDL3::Headers)
target_include_directories(
  imgui_headers
  INTERFACE
    $<BUILD_INTERFACE:${IMGUI_DIR}>
    $<BUILD_INTERFACE:${IMGUI_DIR}/backends>
)
set_target_properties(imgui_headers PROPERTIES PUBLIC_HEADER "${IMGUI_HEADERS}")