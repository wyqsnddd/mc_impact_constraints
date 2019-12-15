#pragma once

#include <math.h>
#include <mc_prediction/mi_qpEstimator.h>
#include <mc_rbdyn/Robots.h>
#include <mc_solver/InequalityConstraint.h>

#include <RBDyn/CoM.h>

#include "ConstraintUtils.h"

namespace mc_impact
{

template<typename Point>
struct DCMWithImpulse : public mc_solver::InequalityConstraintRobot
{

  DCMWithImpulse(mi_qpEstimator & predictor,
                 const mc_rbdyn::Robot & realRobot,
                 const ImpactAwareConstraintParams<Point> & params);

  inline int maxInEq() const override
  {
    return static_cast<int>(getParams().dcmAreaVertexSet.size());
  }

  inline std::string nameInEq() const override
  {
    return "DCMWithImpulse";
  }

  inline const Eigen::MatrixXd & A() const override
  {
    return A_;
  }
  inline const Eigen::VectorXd & bInEq() const override
  {
    return b_;
  }

  void compute() override;

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

  Eigen::Vector3d ComVel_;

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
  inline const std::vector<Point> & getVertices()
  {
    return getParams().dcmAreaVertexSet;
  }
  inline const ImpactAwareConstraintParams<Point> & getParams() const
  {
    return params_;
  }

private:
  // Predictor
  mi_qpEstimator & predictor_;

  const mc_rbdyn::Robot & realRobot_;
  const ImpactAwareConstraintParams<Point> & params_;

  // Alpha vector
  Eigen::VectorXd alpha_;

  ///< Omega
  double omega_ = 1.0;

  // bool debug_;

  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  Eigen::Vector3d predictedComVelJump_;
  std::shared_ptr<rbd::CoMJacobian> comJacobianPtr_;
  // Eigen::MatrixXd comJacobian_;
}; // end of dcmWithImpuplse

} // namespace mc_impact
