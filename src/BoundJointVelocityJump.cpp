#include "BoundJointVelocityJump.h"

namespace mc_impact
{

BoundJointVelocityJump::BoundJointVelocityJump(std::shared_ptr<mi_qpEstimator> predictorPtr, double dt, double multi, bool debug)
: BoundJointVelocityJump(predictorPtr,
                         dt,
                         multi * rbd::dofToVector(predictorPtr->getSimRobot().mb(), predictorPtr->getSimRobot().vl()),
                         multi * rbd::dofToVector(predictorPtr->getSimRobot().mb(), predictorPtr->getSimRobot().vu()),
                         debug)
{
}

BoundJointVelocityJump::BoundJointVelocityJump(
		std::shared_ptr<mi_qpEstimator> predictorPtr,
                                               double dt,
                                               const Eigen::VectorXd & LBound,
                                               const Eigen::VectorXd & UBound,
                                               bool debug)
: mc_solver::GenInequalityConstraintRobot(predictorPtr->getSimRobot().robotIndex()), predictorPtr_(predictorPtr), dt_(dt),
  alpha_L_(LBound), alpha_U_(UBound), debug_(debug)
{
  alpha_.resize(alpha_L_.size());
  if(getPredictor()->getSimRobot().mb().joint(0).dof() == 6)
  {
    startIndex_ = 6;
  }
  int nDof = getPredictor()->getSimRobot().mb().nrDof();
  int realDof = nDof - startIndex_;
  alpha_L_ = alpha_L_.tail(realDof).eval();
  alpha_U_ = alpha_U_.tail(realDof).eval();
  A_.resize(realDof, nDof);
  L_.resize(realDof);
  U_.resize(realDof);

  diff_lower_.resize(realDof);
  diff_upper_.resize(realDof);
  test_delta_vel_.resize(realDof);
}

int BoundJointVelocityJump::maxGenInEq() const
{
  return static_cast<int>(alpha_L_.size());
}

void BoundJointVelocityJump::compute()
{
  const auto & robot = getPredictor()->getSimRobot();
  const auto & J_delta = getPredictor()->getJacobianDeltaAlpha();

  A_ = J_delta.block(startIndex_, 0, J_delta.rows() - startIndex_, J_delta.cols()) * dt_;
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  L_ = alpha_L_ - alpha_.tail(alpha_L_.size())
       - J_delta.block(startIndex_, 0, J_delta.rows() - startIndex_, J_delta.cols()) * alpha_;
  U_ = L_ - alpha_L_ + alpha_U_;

  if(debug_)
  {

    Eigen::VectorXd alphaD = (rbd::dofToVector(getPredictor()->getSimRobot().mb(), getPredictor()->getSimRobot().mbc().alphaD));

    diff_upper_ = U_ - A_ * alphaD;
    diff_lower_ = U_ - diff_upper_ - L_;

    test_delta_vel_ = J_delta * (alphaD * dt_ + alpha_);
  }
}

} // namespace mc_impact
