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
  list(APPEND constraints_SRC 
	  ${CMAKE_SOURCE_DIR}/include/${CSTR}.h 
	  ${CMAKE_SOURCE_DIR}/src/${CSTR}.cpp)
endforeach()

# add_definitions(-fPIC)
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# set(CMAKE_CXX_FLAGS "-Wall -fPIC -o2")

add_library(mc_impact_constraints SHARED ${constraints_SRC})
target_include_directories(mc_impact_constraints PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

target_link_libraries(mc_impact_constraints mc_rtc::mc_control)
target_link_libraries(mc_impact_constraints mc_impact_predictor::mc_impact_predictor)
target_link_libraries(mc_impact_constraints mc_dynamicStability::mc_dynamicStability)
# target_link_libraries(mc_impact_constraints qhullcpp qhullstatic_r)

install(
  TARGETS mc_impact_constraints 
  EXPORT "${TARGETS_EXPORT_NAME}"
  DESTINATION ${MC_RTC_LIBDIR})

install(
	DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION include/mc_impact_constraints
	FILES_MATCHING PATTERN "*.h*")

