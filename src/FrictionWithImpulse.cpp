#include "FrictionWithImpulse.h"

namespace mc_impact
{
FrictionWithImpulse::FrictionWithImpulse(
    const std::string & mcContactName,
    mi_qpEstimator & predictor,
    const ImpactAwareConstraintParams<Eigen::Vector2d> & impactAwareConstraintParams)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), mcContactName_(mcContactName),
  predictor_(predictor), constraintParams_(impactAwareConstraintParams)
{

  // Eigen::Vector3d normal = contact.X_0_r2s(solver.robots()).rotation().row(2).transpose();
  // multiplier_ = (Eigen::MatrixXd::Identity(3, 3) - (1 + mu) * normal * normal.transpose());
  multiplier_.resize(2, 3);
  multiplier_.setZero();
  multiplier_(0, 0) = 1.0;
  multiplier_(0, 2) = -getParams().frictionCoe;
  multiplier_(1, 1) = 1.0;
  multiplier_(1, 2) = -getParams().frictionCoe;

  int nDof = predictor_.getSimRobot().mb().nrDof();

  alpha_.resize(nDof);
  A_.resize(2, nDof);
  b_.resize(2);
}
void FrictionWithImpulse::compute()
{

  const auto & robot = predictor_.getSimRobot();
  const auto & J_deltaF = predictor_.getJacobianDeltaF(getParams().bodyName);
  // std::cout<<"size of J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;
  // std::cout<<"size of reduced J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;

  A_ = (getConstraintParams().dt / getConstraintParams().impactDuration) * multiplier_ * J_deltaF;

  // std::cout<<"size of A_: "<<A_.rows()<<", "<<A_.cols()<<std::endl;
  rbd::paramToVector(robot.mbc().alpha, alpha_);

  // std::cout<<"size of alpha_"<<alpha_.rows()<<std::endl;
  // b_ = -multiplier_ * (predictor_.getSimRobot().forceSensor(sName_).wrench().force() + J_deltaF * alpha_ /
  // impact_dt_);
  b_ = -multiplier_
       * (predictor_.getSimRobot().bodyWrench(getParams().bodyName).force()
          + J_deltaF * alpha_ / getConstraintParams().impactDuration);
}

} // namespace mc_impact
