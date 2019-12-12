#pragma once

#include <Eigen/Dense>
#include <McDynamicStability/McContact.h>
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
  ///< Update the Mc-ZMP-area
  bool updateMcZMPArea = true;
  ///<
  double lowerSlope = 0.01;
  ///<
  double upperSlope = 100.0;
  ///< Debug mode ore not.
  bool debug = false;

  const std::vector<McContactParams> contacts;
  ///< Vertices of the multi-contact DCM area.
  std::vector<Point> dcmAreaVertexSet;
  ///< Vertices of the multi-contact ZMP area.
  std::vector<Point> zmpAreaVertexSet;
};

template struct ImpactAwareConstraintParams<Eigen::Vector3d>;
template struct ImpactAwareConstraintParams<Eigen::Vector2d>;

struct ZMPArea
{
  double min_x;
  double max_x;
  double min_y;
  double max_y;
};

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

} // namespace mc_impact
