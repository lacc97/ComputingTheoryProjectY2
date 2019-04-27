cmake_minimum_required(VERSION 3.12)
project(computingproject)

set(CMAKE_CXX_STANDARD 14)

find_package(Python3    REQUIRED COMPONENTS Development)
find_package(Eigen3     REQUIRED NO_MODULE)
find_package(OpenMP     REQUIRED)

add_executable(computingproject main.cpp utils/timer.cpp utils/types.cpp utils/utils.hpp)
    target_compile_definitions(computingproject PUBLIC WITHOUT_NUMPY)
    target_include_directories(computingproject PRIVATE utils/)
    target_link_libraries(computingproject Python3::Python Eigen3::Eigen OpenMP::OpenMP_CXX)

add_subdirectory(test/)
