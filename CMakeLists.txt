cmake_minimum_required(VERSION 2.8.12)
include(CheckCXXCompilerFlag)

# -march
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

# project
project(py_power_diagram_2d_double)
add_subdirectory(ext/pybind11)
pybind11_add_module(py_power_diagram_2d_double py_power_diagram/py_power_diagram_2d_double.cpp)
target_include_directories(py_power_diagram_2d_double PRIVATE ext/xsimd/include)
set_property(TARGET py_power_diagram_2d_double PROPERTY CXX_STANDARD 14)
