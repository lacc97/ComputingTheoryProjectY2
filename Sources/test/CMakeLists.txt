find_package(PkgConfig REQUIRED)
    pkg_search_module(SDL2 REQUIRED sdl2)

add_executable(double_pendulum double_pendulum.cpp ../utils/types.cpp ../utils/timer.cpp)
    target_include_directories(double_pendulum PRIVATE ../utils/)
    target_link_libraries(double_pendulum Eigen3::Eigen sdl2pp)

set(ENABLE_SDL2_GFX)
add_subdirectory(sdl2/)
