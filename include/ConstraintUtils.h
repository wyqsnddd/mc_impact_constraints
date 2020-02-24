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
public:
  ImpactAwareConstraintParams()
  {
    //contactSetPtr.reset(new McContactSet());
  }

  ~ImpactAwareConstraintParams() {}
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
