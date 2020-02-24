#include "BoundJointTorqueJump.h"

namespace mc_impact
{

BoundJointTorqueJump::BoundJointTorqueJump(
		std::shared_ptr<mi_qpEstimator> predictorPtr,
                                           double dt,
                                           double impact_dt,
                                           double mult,
                                           bool debug)
: BoundJointTorqueJump(predictorPtr,
                       dt,
                       impact_dt,
                       mult * rbd::dofToVector(predictorPtr->getSimRobot().mb(), predictorPtr->getSimRobot().tl()),
                       mult * rbd::dofToVector(predictorPtr->getSimRobot().mb(), predictorPtr->getSimRobot().tu()),
                       debug)
{
}

BoundJointTorqueJump::BoundJointTorqueJump(std::shared_ptr<mi_qpEstimator> predictorPtr,
                                           double dt,
                                           double impact_dt,
                                           const Eigen::VectorXd & LBound,
                                           const Eigen::VectorXd & UBound,
                                           bool debug)
: mc_solver::GenInequalityConstraintRobot(predictorPtr->getSimRobot().robotIndex()), predictorPtr_(predictorPtr), dt_(dt),
  impact_dt_(impact_dt), tau_L_(LBound), tau_U_(UBound), debug_(debug)
{
  // std::cout<<"The lower impulsive torque bound is: "<<std::endl<<LBound.transpose()<<std::endl;
  // std::cout<<"The upper impulsive torque bound is: "<<std::endl<<UBound.transpose()<<std::endl;
  alpha_.resize(tau_L_.size());
  if(getPredictor()->getSimRobot().mb().joint(0).dof() == 6)
  {
    startIndex_ = 6;
  }
  int nDof = getPredictor()->getSimRobot().mb().nrDof();
  int realDof = nDof - startIndex_;
  tau_L_ = tau_L_.tail(realDof).eval();
  tau_U_ = tau_U_.tail(realDof).eval();
  A_.resize(realDof, nDof);
  L_.resize(realDof);
  U_.resize(realDof);
  test_delta_torque_.resize(realDof);
  difference_lower_.resize(realDof);
  difference_upper_.resize(realDof);
}

int BoundJointTorqueJump::maxGenInEq() const
{
  return static_cast<int>(tau_L_.size());
}

void BoundJointTorqueJump::compute()
{
  const auto & robot = getPredictor()->getSimRobot();
  const auto & J_deltatau = getPredictor()->getJacobianDeltaTau();
  A_ = (dt_ / impact_dt_) * J_deltatau.block(startIndex_, 0, J_deltatau.rows() - startIndex_, J_deltatau.cols());
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  L_ = tau_L_
       - J_deltatau.block(startIndex_, 0, J_deltatau.rows() - startIndex_, J_deltatau.cols()) * alpha_ / impact_dt_;
  U_ = L_ - tau_L_ + tau_U_;

  if(debug_)
  {

    Eigen::VectorXd alphaD = (rbd::dofToVector(getPredictor()->getSimRobot().mb(), getPredictor()->getSimRobot().mbc().alphaD));

    difference_upper_ = U_ - A_ * alphaD;
    difference_lower_ = U_ - difference_upper_ - L_;

    test_delta_torque_ = (1.0 / impact_dt_) * J_deltatau * (alpha_ + alphaD * dt_);
  }
}

} // namespace mc_impact
