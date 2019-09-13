#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_rbdyn/Robots.h>
#include <mc_solver/InequalityConstraint.h>

#include <RBDyn/CoM.h>

#include "constraintUtils.h"
#include <math.h>

namespace mc_impact
{

template<typename supportContact, typename Point>
struct dcmWithImpulse : public mc_solver::InequalityConstraint
{

  dcmWithImpulse(const mc_rbdyn::Robot & realRobot,
                 mi_qpEstimator & predictor,
                 const std::vector<supportContact> & supports,
                 double dt,
                 double impact_dt,
                 const std::vector<Point> & vertexSet,
                 double lowerSlope = 0.01,
                 double upperSlope = 100,
                 bool debug = false);

  inline int maxInEq() const override
  {
    return static_cast<int>(iniVertexSet_.size());
  }

  inline std::string nameInEq() const override
  {
    return "dcmWithImpulse";
  }

  inline const Eigen::MatrixXd & A() const override
  {
    return A_;
  }
  inline const Eigen::VectorXd & bInEq() const override
  {
    return b_;
  }

  void computeAb() override;

  inline const Eigen::MatrixXd & getA()
  {
    return A_;
  }
  inline const Eigen::VectorXd & getb()
  {
    return b_;
  }

  // Debugging:
  Point centeroid_;
  Eigen::VectorXd slopeVec_ = Eigen::Vector3d::Zero();

  bool zmpTest_ = false;
  Eigen::MatrixXd G_dcm_;
  Eigen::VectorXd h_dcm_;

  Eigen::Vector2d dcm_;
  Eigen::Vector2d predicted_dcm_;
  Eigen::Vector2d predicted_dcm_jump_;
  Eigen::VectorXd difference_;

  Eigen::MatrixXd A_dcm_;
  bool pointInsideSupportPolygon(const Point & input);

  inline void calcOmega(const double & c_z)
  {
    assert(c_z > 0.0);
    omega_ = sqrt(9.8 / c_z);
  }
  inline const double & getOmega()
  {
    return omega_;
  }
  inline const Eigen::Vector3d & getPredictedComVelJump()
  {
    return predictedComVelJump_;
  }

private:
  const mc_rbdyn::Robot & realRobot_;
  // Predictor
  mi_qpEstimator & predictor_;
  // Timestep
  double dt_;
  // Impact duration
  double impact_dt_;
  // Alpha vector
  Eigen::VectorXd alpha_;

  const std::vector<Point> iniVertexSet_;
  double omega_ = 1.0;

  bool debug_;

  std::vector<supportContact> supports_;
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  Eigen::Vector3d predictedComVelJump_;
  std::shared_ptr<rbd::CoMJacobian> comJacobianPtr_;
  // Eigen::MatrixXd comJacobian_;
}; // end of dcmWithImpuplse

} // namespace mc_impact
