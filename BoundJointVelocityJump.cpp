#include "BoundJointVelocityJump.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

BoundJointVelocityJump::BoundJointVelocityJump(mi_impactPredictor & predictor,
                                               double dt)
: BoundJointVelocityJump(predictor, dt,
                         rbd::dofToVector(predictor.getRobot().mb(), predictor.getRobot().vl()),
                         rbd::dofToVector(predictor.getRobot().mb(), predictor.getRobot().vu()))
{
}

BoundJointVelocityJump::BoundJointVelocityJump(mi_impactPredictor & predictor,
                                               double dt,
                                               const Eigen::VectorXd & LBound,
                                               const Eigen::VectorXd & UBound)
: mc_solver::GenInequalityConstraint(predictor.getRobot().robotIndex()),
  predictor_(predictor),
  dt_(dt), alpha_L_(LBound), alpha_U_(UBound)
{
  alpha_.resize(alpha_L_.size());
}

int BoundJointVelocityJump::maxGenInEq() const
{
  return predictor_.getRobot().mb().nrDof();
}

void BoundJointVelocityJump::computeALU()
{
  const auto & robot = predictor_.getRobot();
  const auto & J_delta = predictor_.getJacobianDeltaAlpha();
  A_ = J_delta / dt_;
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  L_ = alpha_L_ - alpha_ - J_delta * alpha_;
  U_ = L_ - alpha_L_ + alpha_U_;
}

}
