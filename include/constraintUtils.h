# pragma once

#include <Eigen/Dense>
#include <iostream>

namespace mc_impact
{

struct ZMPArea
{
  double min_x;
  double max_x;
  double min_y;
  double max_y;
};
struct zmpSupportContact
{
  std::string bodyName;
  std::string sensorName;
};



template <typename Point>
void pointsToInequalityMatrix(const std::vector<Point> & inputPoints, Eigen::MatrixXd & G, Eigen::VectorXd & h, Point & centeroid, Eigen::VectorXd & slopeVec, double miniSlope = 0.01, double maxSlope = 100);

template <typename Point>
void pointsToInequalityMatrix(const std::vector<Point> & inputPoints, Eigen::MatrixXd & G, Eigen::VectorXd & h, double miniSlope = 0.01, double maxSlope = 100);


}// end of namespace 
