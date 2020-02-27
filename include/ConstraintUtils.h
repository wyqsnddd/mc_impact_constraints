#pragma once

#include <Eigen/Dense>
#include <McDynamicStability/McContact.h>
#include <McDynamicStability/Utils.h>
#include <iostream>

namespace mc_impact
{

template<typename Point>
struct ImpactAwareConstraintParams
{
  ///< Sampling Period
  double dt = 0.05;
  ///< Impact duration
  double impactDuration = 0.05;
  ///< Consider multiple contacts (not only the two feet contacts)
  bool multiContactCase = true;

  ///< Update the Mc-ZMP-area for the ZMP constraint
  bool updateMcZMPArea = true;
  ///< Update the Mc-DCM-area for the DCM constraint
  bool updateMcDCMArea = true;

  ///< Debug mode ore not.
  bool debug = false;

  ///<   ProjectionSlop
  double lowerSlope = 0.01;
  ///<   ProjectionSlop
  double upperSlope = 100.0;

  ///< Update the Mc-ZMP-area
  bool constrainingZMP = false;
  bool constrainingDCM = false;
  bool constrainingCOMAcc = false;

  //std::shared_ptr<McContactSet> contactSetPtr;
  // const std::vector<McContactParams> contacts;
  ///< Vertices of the multi-contact DCM area.
  std::vector<Point> dcmAreaVertexSet;
  ///< Vertices of the multi-contact ZMP area.
  std::vector<Point> zmpAreaVertexSet;

  ///< Specify the projection of McZMP and McCom areas.
  McProjectionParams mcProjectionParams;

  bool enableImpactAwareFloatingBaseConstraint = false;
  bool enableImpactAwareJointVelocityConstraint = false;
  bool enableImpactAwareJointTorqueConstraint = false;
  bool enableImpactAwareCOPConstraint = false;
  bool enableImpactAwareFrictionConeConstraint = false;
};

struct ImpactAwareState
{
  Eigen::Vector3d current = Eigen::Vector3d::Zero();
  Eigen::Vector3d stateJump = Eigen::Vector3d::Zero();
  Eigen::Vector3d oneStepPreview = Eigen::Vector3d::Zero();
};
/*
struct ComparisonState
{
  Eigen::Vector3d realRobot;
  Eigen::Vector3d simRobot;
};
*/
struct FloatingBaseStates
{
  // These values are all taken from the realrobot or the estimator
  ImpactAwareState mcZMP;
  ImpactAwareState bipedalZMP;
  ImpactAwareState ComVel;
  ImpactAwareState DCM;

  Eigen::Vector3d Com;
  Eigen::Vector3d ComAcc;
  // Eigen::Vector3d ComVelSimRobot;
};

template struct ImpactAwareConstraintParams<Eigen::Vector3d>;
template struct ImpactAwareConstraintParams<Eigen::Vector2d>;

/*
template<typename Point>
void pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                              Eigen::MatrixXd & G,
                              Eigen::VectorXd & h,
                              Point & centeroid,
                              Eigen::VectorXd & slopeVec,
                              double miniSlope = 0.01,
                              double maxSlope = 100);
template<typename Point>
void pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                              Eigen::MatrixXd & G,
                              Eigen::VectorXd & h,
                              double miniSlope = 0.01,
                              double maxSlope = 100);

            */
} // namespace mc_impact

/*
// Reader of McContactParams
namespace mc_rtc
{
template<typename Point>
struct ConfigurationLoader<mc_impact::ImpactAwareConstraintParams<Point> >
{
  static mc_impact::ImpactAwareConstraintParams<Point> load(const mc_rtc::Configuration & config)
  {
    mc_impact::ImpactAwareConstraintParams<Point> impactAwareConstraintParams;
    impactAwareConstraintParams.dt = config("impact")("timeStep");
  
    impactAwareConstraintParams.impactDuration = config("impact")("impactDuration");
  
    impactAwareConstraintParams.updateMcZMPArea =
        static_cast<bool>(config("impact")("constraints")("floatingBaseConstraint")("zmpArea")("updateMcZMPArea"));
    impactAwareConstraintParams.updateMcDCMArea =
        static_cast<bool>(config("impact")("constraints")("floatingBaseConstraint")("dcmArea")("updateMcDCMArea"));
  
     impactAwareConstraintParams.debug =  static_cast<bool>(config("impact")("constraints")("debug"));
    impactAwareConstraintParams.dcmAreaVertexSet = config("impact")("constraints")("initialDCMArea");
    impactAwareConstraintParams.zmpAreaVertexSet = config("impact")("constraints")("initialZMPArea");
  
    impactAwareConstraintParams.enableImpactAwareFloatingBaseConstraint =
        static_cast<bool>(config("impact")("constraints")("floatingBaseConstraint")("enabled"));
  
    impactAwareConstraintParams.enableImpactAwareJointVelocityConstraint =
        static_cast<bool>(config("impact")("constraints")("jointVelocity")("on"));
    impactAwareConstraintParams.enableImpactAwareJointTorqueConstraint =
        static_cast<bool>(config("impact")("constraints")("jointTorque")("on"));



    impactAwareConstraintParams.mcProjectionParams.iterationLimit =
          config("impact")("constraints")("floatingBaseConstraint")("mcProjectionParams")("projectionIterationLimit");
      impactAwareConstraintParams.mcProjectionParams.convergeThreshold = config("impact")("constraints")(
          "floatingBaseConstraint")("mcProjectionParams")("projectionConvergenceThreshold");
      impactAwareConstraintParams.mcProjectionParams.projectionRadius =
          config("impact")("constraints")("floatingBaseConstraint")("mcProjectionParams")("projectionRadius");
      impactAwareConstraintParams.mcProjectionParams.useLIPMAssumptions = static_cast<bool>(
          config("impact")("constraints")("floatingBaseConstraint")("mcProjectionParams")("useLIPMAssumptions"));
      impactAwareConstraintParams.mcProjectionParams.useSpatialVectorAlgebra = static_cast<bool>(
          config("impact")("constraints")("floatingBaseConstraint")("mcProjectionParams")("useSpatialVectorAlgebra"));
      impactAwareConstraintParams.mcProjectionParams.debug = static_cast<bool>(
          config("impact")("constraints")("floatingBaseConstraint")("debug"));

      impactAwareConstraintParams.constrainingDCM =
          static_cast<bool>(config("impact")("constraints")("floatingBaseConstraint")("dcmArea")("dcmConstraint"));
      impactAwareConstraintParams.constrainingZMP =
          static_cast<bool>(config("impact")("constraints")("floatingBaseConstraint")("zmpArea")("zmpConstraint"));
      impactAwareConstraintParams.constrainingCOMAcc =
          static_cast<bool>(config("impact")("constraints")("floatingBaseConstraint")("comArea")("comAccConstraint"));


    return impactAwareConstraintParams;

  }

  //template<typename Point>
  static mc_rtc::Configuration save(const mc_impact::ImpactAwareConstraintParams<Point> & data)
  {
    mc_rtc::Configuration config;

    config.add("timeStep", data.dt);
    config.add("impactDuration", data.impactDuration);
    config.add("updateMcZMPArea", data.updateMcZMPArea);
    config.add("updateMcDCMArea", data.updateMcDCMArea);
    config.add("debug", data.debug);
    config.add("dcmAreaVertexSet", data.dcmAreaVertexSet);
    config.add("zmpAreaVertexSet", data.zmpAreaVertexSet);

    // Constraints configuration
    config.add("enableImpactAwareFloatingBaseConstraint", data.enableImpactAwareFloatingBaseConstraint);
    config.add("enableImpactAwareJointVelocityConstraint", data.enableImpactAwareJointVelocityConstraint);
    config.add("enableImpactAwareJointTorqueConstraint", data.enableImpactAwareJointTorqueConstraint);


    // Configuration of the FloatingBaseConstraint
    config.add("constrainingDCM", data.constrainingDCM);
    config.add("constrainingZMP", data.constrainingZMP);
    config.add("constrainingCOMAcc", data.constrainingCOMAcc);


    // Projection parameters 
    config.add("iterationLimit", data.mcProjectionParams.iterationLimit);
    config.add("convergeThreshold", data.mcProjectionParams.convergeThreshold);
    config.add("projectionRadius", data.mcProjectionParams.projectionRadius);
    config.add("useLIPMAssumptions", data.mcProjectionParams.useLIPMAssumptions);
    config.add("useSpatialVectorAlgebra", data.mcProjectionParams.useSpatialVectorAlgebra);
    config.add("debug", data.mcProjectionParams.debug);
    
    return config;
  }
};
{
struct ConfigurationLoader<mc_impact::ImpactAwareConstraintParams<Eigen::Vector2d> >;
struct ConfigurationLoader<mc_impact::ImpactAwareConstraintParams<Eigen::Vector3d> >;

} // namespace mc_rtc 

*/



