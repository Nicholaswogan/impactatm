cmake_minimum_required(VERSION "3.14")
project(impactatm C)
include(cmake/CPM.cmake)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

if (NOT SKBUILD)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
endif()

find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(Cython REQUIRED)
add_subdirectory(impactatm)

if (SKBUILD)
  install(TARGETS _cvode DESTINATION impactatm)
else()
  install(TARGETS _cvode DESTINATION impactatm)
endif()
