cmake_minimum_required(VERSION 3.1)

include(cmake/base.cmake)
include(cmake/boost.cmake)

set(PROJECT_NAME mc_impact_constraints)

set(PROJECT_DESCRIPTION "The QP impact constraints")
set(PROJECT_URL  https://gite.lirmm.fr/jrl-collaborators/mc_impact_constraints )
set(PROJECT_VERSION 1.0.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)


project(${PROJECT_NAME} VERSION ${PROJECT_VERSION})
setup_project()

find_package(mc_rtc REQUIRED)
find_package(mc_impact_predictor REQUIRED)
find_package(mc_dynamicStability REQUIRED)
find_package(Stabiliplus REQUIRED)

include_directories(include)

add_subdirectory(src)

install(
  FILES "${CMAKE_SOURCE_DIR}/mc_impact_constraintsConfig.cmake"
  DESTINATION lib/cmake/mc_impact_constraints
)

setup_project_finalize()
