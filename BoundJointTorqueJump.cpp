#include "BoundJointTorqueJump.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

BoundJointTorqueJump::BoundJointTorqueJump(mi_impactPredictor & predictor, double dt, double impact_dt, double mult)
: BoundJointTorqueJump(predictor,
                       dt,
                       impact_dt,
                       mult * rbd::dofToVector(predictor.getRobot().mb(), predictor.getRobot().tl()),
                       mult * rbd::dofToVector(predictor.getRobot().mb(), predictor.getRobot().tu()))
{
}

BoundJointTorqueJump::BoundJointTorqueJump(mi_impactPredictor & predictor,
                                           double dt,
                                           double impact_dt,
                                           const Eigen::VectorXd & LBound,
                                           const Eigen::VectorXd & UBound)
: mc_solver::GenInequalityConstraint(predictor.getRobot().robotIndex()), predictor_(predictor), dt_(dt),
  impact_dt_(impact_dt), tau_L_(LBound), tau_U_(UBound)
{
  alpha_.resize(tau_L_.size());
  if(predictor_.getRobot().mb().joint(0).dof() == 6)
  {
    startIndex_ = 6;
  }
  int nDof = predictor_.getRobot().mb().nrDof();
  tau_L_ = tau_L_.tail(nDof - startIndex_).eval();
  tau_U_ = tau_U_.tail(nDof - startIndex_).eval();
  A_.resize(nDof - startIndex_, nDof);
  L_.resize(nDof - startIndex_);
  U_.resize(nDof - startIndex_);
}

int BoundJointTorqueJump::maxGenInEq() const
{
  return tau_L_.size();
}

void BoundJointTorqueJump::computeALU()
{
  const auto & robot = predictor_.getRobot();
  const auto & J_deltatau = predictor_.getJacobianDeltaTau();
  A_ = (dt_ / impact_dt_) * J_deltatau.block(startIndex_, 0, J_deltatau.rows() - startIndex_, J_deltatau.cols());
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  L_ = tau_L_
       - J_deltatau.block(startIndex_, 0, J_deltatau.rows() - startIndex_, J_deltatau.cols()) * alpha_ / impact_dt_;
  U_ = L_ - tau_L_ + tau_U_;
}

} // namespace mc_impact
