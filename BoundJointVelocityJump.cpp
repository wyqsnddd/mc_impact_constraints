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
  if(predictor_.getRobot().mb().joint(0).dof() == 6)
  {
    startIndex_ = 6;
  }
  alpha_ = alpha_.tail(alpha_.size() - startIndex_);
  alpha_L_ = alpha_L_.tail(alpha_.size());
  alpha_U_ = alpha_U_.tail(alpha_.size());
}

int BoundJointVelocityJump::maxGenInEq() const
{
  return alpha_.size();
}

void BoundJointVelocityJump::computeALU()
{
  const auto & robot = predictor_.getRobot();
  const auto & J_delta = predictor_.getJacobianDeltaAlpha();
  A_ = J_delta.block(startIndex_, 0, J_delta.rows() - startIndex_, J_delta.cols()) / dt_;
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  alpha_ = alpha_.tail(alpha_L_.size());
  L_ = alpha_L_ - alpha_ - J_delta.block(startIndex_, 0, J_delta.rows() - startIndex_, J_delta.cols()) * alpha_;
  U_ = L_ - alpha_L_ + alpha_U_;
}

}
