set(constraints
  BoundJointVelocityJump
  BoundJointTorqueJump
  ContactWrenchMatrixToLambdaMatrix
  COPInsideContactAreaWithImpulse
  PositiveContactForceWithImpulse
  ZMPInsideSupportPolygonWithImpulse
)

set(constraints_SRC)
foreach(CSTR ${constraints})
  list(APPEND constraints_SRC ${CSTR}.h ${CSTR}.cpp)
endforeach()

add_library(mc_impact_constraints SHARED ${constraints_SRC})
target_include_directories(mc_impact_constraints PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
pkg_config_use_dependency(mc_impact_constraints mc_rtc)
install(TARGETS mc_impact_constraints DESTINATION ${MC_RTC_LIBDIR})
