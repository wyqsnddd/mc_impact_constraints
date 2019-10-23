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

include_directories(include)

set(constraints
  BoundJointVelocityJump
  BoundJointTorqueJump
  #ContactWrenchMatrixToLambdaMatrix
  frictionWithImpulse
  copWithImpulse
  zmpWithImpulse
  dcmWithImpulse
  #supportPolygon
  constraintUtils
)
set(constraints_SRC)
foreach(CSTR ${constraints})
  list(APPEND constraints_SRC include/${CSTR}.h src/${CSTR}.cpp)
endforeach()

# add_definitions(-fPIC)
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# set(CMAKE_CXX_FLAGS "-Wall -fPIC -o2")

add_library(mc_impact_constraints SHARED ${constraints_SRC})
target_include_directories(mc_impact_constraints PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(mc_impact_constraints mc_rtc::mc_control)
target_link_libraries(mc_impact_constraints mc_impact_predictor::mc_impact_predictor)
# target_link_libraries(mc_impact_constraints qhullcpp qhullstatic_r)

install(
  TARGETS mc_impact_constraints 
  EXPORT "${TARGETS_EXPORT_NAME}"
  DESTINATION ${MC_RTC_LIBDIR})

install(
	DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION include/mc_impact_constraints
	FILES_MATCHING PATTERN "*.h*")

setup_project_finalize()
setup_project_package_finalize()
