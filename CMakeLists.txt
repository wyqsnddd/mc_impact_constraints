cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME mc_impact_constraints)
set(PROJECT_DESCRIPTION "The QP impact constraints")
set(PROJECT_URL https://gite.lirmm.fr/jrl-collaborators/mc_impact_constraints)
set(PROJECT_VERSION 1.0.0)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(PROJECT_DEBUG_POSTFIX "_d")
set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)

include(cmake/base.cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})

add_project_dependency(mc_impact_predictor REQUIRED)
add_project_dependency(mc_dynamicStability REQUIRED)

add_subdirectory(src)
