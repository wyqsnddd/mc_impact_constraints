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
pkg_config_use_dependency(mc_impact_constraints mc_rtc )
target_link_libraries(mc_impact_constraints mc_impact_predictor)


pkg_config_append_libs(mc_impact_constraints)

install(TARGETS mc_impact_constraints DESTINATION ${MC_RTC_LIBDIR})

install(
	DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION include/mc_impact_constraints
	FILES_MATCHING PATTERN "*.h*")




