#include "BoundJointVelocityJump.h"

namespace mc_impact
{

BoundJointVelocityJump::BoundJointVelocityJump(mi_qpEstimator & predictor, double dt, bool debug)
: BoundJointVelocityJump(predictor,
                         dt,
                         rbd::dofToVector(predictor.getSimRobot().mb(), predictor.getSimRobot().vl()),
                         rbd::dofToVector(predictor.getSimRobot().mb(), predictor.getSimRobot().vu()),
                         debug)
{
}

BoundJointVelocityJump::BoundJointVelocityJump(mi_qpEstimator & predictor,
                                               double dt,
                                               const Eigen::VectorXd & LBound,
                                               const Eigen::VectorXd & UBound,
                                               bool debug)
: mc_solver::GenInequalityConstraint(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt),
  alpha_L_(LBound), alpha_U_(UBound), debug_(debug)
{
  alpha_.resize(alpha_L_.size());
  if(predictor_.getSimRobot().mb().joint(0).dof() == 6)
  {
    startIndex_ = 6;
  }
  int nDof = predictor_.getSimRobot().mb().nrDof();
  alpha_L_ = alpha_L_.tail(nDof - startIndex_).eval();
  alpha_U_ = alpha_U_.tail(nDof - startIndex_).eval();
  A_.resize(nDof - startIndex_, nDof);
  L_.resize(nDof - startIndex_);
  U_.resize(nDof - startIndex_);

  diff_lower_.resize(nDof - startIndex_);
  diff_upper_.resize(nDof - startIndex_);
}

int BoundJointVelocityJump::maxGenInEq() const
{
  return alpha_L_.size();
}

void BoundJointVelocityJump::computeALU()
{
  const auto & robot = predictor_.getSimRobot();
  const auto & J_delta = predictor_.getJacobianDeltaAlpha();
  A_ = J_delta.block(startIndex_, 0, J_delta.rows() - startIndex_, J_delta.cols()) * dt_;
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  L_ = alpha_L_ - alpha_.tail(alpha_L_.size())
       - J_delta.block(startIndex_, 0, J_delta.rows() - startIndex_, J_delta.cols()) * alpha_;
  U_ = L_ - alpha_L_ + alpha_U_;

  if(debug_)
  {

    Eigen::VectorXd alphaD = (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD));

    diff_upper_ = U_ - A_ * alphaD;
    diff_lower_ = U_ - diff_upper_ - L_;
  }
}

} // namespace mc_impact
