find_dependency(mc_impact_predictor REQUIRED)
find_dependency(mc_dynamicStability REQUIRED)
find_dependency(Stabiliplus REQUIRED)

if(NOT TARGET mc_impact_constraints::mc_impact_constraints)
  include("${CMAKE_CURRENT_LIST_DIR}/mc_impact_constraintsTargets.cmake")
endif()
