#include "constraintUtils.h"

template<typename T>
T sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

void clampSlope(double & slope, const double & mini, const double & max)
{
  double sign = sgn<double>(slope);

  if(fabs(slope) <= mini)
  {
    slope = sign * mini;
  }
  else if(fabs(slope) >= max)
  {
    slope = sign * max;
  }
}

template<typename Point>
void mc_impact::pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                                         Eigen::MatrixXd & G,
                                         Eigen::VectorXd & h,
                                         Point & center,
                                         Eigen::VectorXd & slopeVec,
                                         double miniSlope,
                                         double maxSlope)
{

  int vertexNumber = static_cast<int>(inputPoints.size());
  int dim = static_cast<int>(inputPoints[0].size());

  // Point center;
  center[0] = 0;
  center[1] = 0;

  G.resize(vertexNumber, dim);
  h.resize(vertexNumber);
  G.setOnes();
  h.setOnes();

  for(auto & p : inputPoints)
  {
    center[0] += p[0];
    center[1] += p[1];
  }
  center[0] = center[0] / (double)vertexNumber;
  center[1] = center[1] / (double)vertexNumber;
  std::cout << "The initial center is: " << center.transpose() << std::endl;

  int vNumber = 0;
  slopeVec = Eigen::VectorXd::Zero(vertexNumber);

  for(auto idx = inputPoints.begin(); idx != inputPoints.end(); idx++, vNumber++)
  {

    Point point_one = *idx;
    Point point_two;

    if((idx + 1) == inputPoints.end())
    {
      point_two = *(inputPoints.begin());
    }
    else
    {
      point_two = *(idx + 1);
    }

    Point difference = point_two - point_one;
    //difference.normalize();
    double slope = difference[1] / difference[0];

    clampSlope(slope, miniSlope, maxSlope);

    G(vNumber, 0) = -slope;
    h(vNumber) = -slope * point_one[0] + point_one[1];

    slopeVec(vNumber) = slope;

    int lineSign = 1;
    /// should remove slope?
    if(!((-slope * center[0] + center[1]) <= h(vNumber)))
    {
      lineSign = -1;
    }
    // Correct the sign with the centeroid point
    // G.block(vNumber, 0, 1, 2) = lineSign * G.block(vNumber, 0, 1, 2);
    // h(vNumber) = lineSign*h(vNumber);

    G.block(vNumber, 0, 1, 2) *= lineSign;
    h(vNumber) *= lineSign;

  } // end of iterating over points

} // end of pointsToInequalityMatrix

template void mc_impact::pointsToInequalityMatrix<Eigen::Vector2d>(const std::vector<Eigen::Vector2d> & inputPoints,
                                                                   Eigen::MatrixXd & G,
                                                                   Eigen::VectorXd & h,
                                                                   Eigen::Vector2d & center,
                                                                   Eigen::VectorXd & slopeVec,
                                                                   double miniSlope,
                                                                   double maxSlope);
template void mc_impact::pointsToInequalityMatrix<Eigen::Vector3d>(const std::vector<Eigen::Vector3d> & inputPoints,
                                                                   Eigen::MatrixXd & G,
                                                                   Eigen::VectorXd & h,
                                                                   Eigen::Vector3d & center,
                                                                   Eigen::VectorXd & slopeVec,
                                                                   double miniSlope,
                                                                   double maxSlope);

template<typename Point>
void mc_impact::pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                                         Eigen::MatrixXd & G,
                                         Eigen::VectorXd & h,
                                         double miniSlope,
                                         double maxSlope)
{

  int vertexNumber = static_cast<int>(inputPoints.size());
  int dim = static_cast<int>(inputPoints[0].size());

  Point center;
  center[0] = 0;
  center[1] = 0;

  G.resize(vertexNumber, dim);
  h.resize(vertexNumber);
  G.setOnes();
  h.setOnes();

  for(auto & p : inputPoints)
  {
    center[0] += p[0];
    center[1] += p[1];
  }
  center[0] = center[0] / (double)vertexNumber;
  center[1] = center[1] / (double)vertexNumber;

  int vNumber = 0;

  for(auto idx = inputPoints.begin(); idx != inputPoints.end(); idx++, vNumber++)
  {

    Point point_one = *idx;
    Point point_two;

    if((idx + 1) == inputPoints.end())
    {
      point_two = *(inputPoints.begin());
    }
    else
    {
      point_two = *(idx + 1);
    }

    Point difference = point_two - point_one;
    difference.normalize();
    double slope = difference[1] / difference[0];

    G(vNumber, 0) = -slope;
    h(vNumber) = -slope * point_one[0] + point_one[1];

    int lineSign = 1;
    /// should remove slope?
    if(!((center[1] - slope * center[0]) <= h(vNumber)))
    {
      lineSign = -1;
    }
    // Correct the sign with the centeroid point
    G.block(vNumber, 0, 1, 2) *= lineSign;
    h(vNumber) *= lineSign;

  } // end of iterating over points

} // end of pointsToInequalityMatrix

template void mc_impact::pointsToInequalityMatrix<Eigen::Vector2d>(const std::vector<Eigen::Vector2d> & inputPoints,
                                                                   Eigen::MatrixXd & G,
                                                                   Eigen::VectorXd & h,
                                                                   double miniSlope,
                                                                   double maxSlope);
template void mc_impact::pointsToInequalityMatrix<Eigen::Vector3d>(const std::vector<Eigen::Vector3d> & inputPoints,
                                                                   Eigen::MatrixXd & G,
                                                                   Eigen::VectorXd & h,
                                                                   double miniSlope,
                                                                   double maxSlope);
