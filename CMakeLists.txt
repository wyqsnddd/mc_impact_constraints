cmake_minimum_required(VERSION 3.1)

include(cmake/base.cmake)
include(cmake/boost.cmake)

set(PROJECT_NAME mc_impact_constraints)

set(PROJECT_DESCRIPTION "The QP impact constraints")
set(PROJECT_URL  https://gite.lirmm.fr/jrl-collaborators/mc_impact_constraints )
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)
setup_project()

find_package(mc_rtc REQUIRED)

include_directories(include)

set(constraints
  BoundJointVelocityJump
  BoundJointTorqueJump
  ContactWrenchMatrixToLambdaMatrix
  ZeroSlippageWithImpulse
  COPInsideContactAreaWithImpulse
  frictionWithImpulse
  copWithImpulse
  zmpWithImpulse
  #ZMPInsideSupportPolygonWithImpulse
)
set(constraints_SRC)
foreach(CSTR ${constraints})
  list(APPEND constraints_SRC include/${CSTR}.h src/${CSTR}.cpp)
endforeach()

add_library(mc_impact_constraints SHARED ${constraints_SRC})
target_include_directories(mc_impact_constraints PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(mc_impact_constraints mc_rtc::mc_control)
target_link_libraries(mc_impact_constraints mc_impact_predictor)

pkg_config_append_libs(mc_impact_constraints)

install(TARGETS mc_impact_constraints DESTINATION ${MC_RTC_LIBDIR})

install(
	DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION include/mc_impact_constraints
	FILES_MATCHING PATTERN "*.h*")

setup_project_finalize()
