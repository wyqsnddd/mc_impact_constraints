#include "FrictionWithImpulse.h"

namespace mc_impact
{

FrictionWithImpulse::FrictionWithImpulse(mi_qpEstimator & predictor,
                                         const std::string & bodyName,
                                         const std::string & sensorName,
                                         //const mc_rbdyn::Contact & contact,
                                         double dt,
                                         double impact_dt,
                                         double mu)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt), impact_dt_(impact_dt)
{

  // Eigen::Vector3d normal = contact.X_0_r2s(solver.robots()).rotation().row(2).transpose();
  // multiplier_ = (Eigen::MatrixXd::Identity(3, 3) - (1 + mu) * normal * normal.transpose());
  multiplier_.resize(2, 3);
  multiplier_.setZero();
  multiplier_(0, 0) = 1.0;
  multiplier_(0, 2) = -mu;
  multiplier_(1, 1) = 1.0;
  multiplier_(1, 2) = -mu;

  bName_ = bodyName;
  sName_ = sensorName;

  int nDof = predictor_.getSimRobot().mb().nrDof();

  alpha_.resize(nDof);
  A_.resize(2, nDof);
  b_.resize(2);
}
void FrictionWithImpulse::compute()
{

  const auto & robot = predictor_.getSimRobot();
  const auto & J_deltaF = predictor_.getJacobianDeltaF(bName_);
  // std::cout<<"size of J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;
  // std::cout<<"size of reduced J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;

  A_ = (dt_ / impact_dt_) * multiplier_ * J_deltaF;

  // std::cout<<"size of A_: "<<A_.rows()<<", "<<A_.cols()<<std::endl;
  rbd::paramToVector(robot.mbc().alpha, alpha_);

  // std::cout<<"size of alpha_"<<alpha_.rows()<<std::endl;
  // b_ = -multiplier_ * (predictor_.getSimRobot().forceSensor(sName_).wrench().force() + J_deltaF * alpha_ /
  // impact_dt_);
  b_ = -multiplier_ * (predictor_.getSimRobot().bodyWrench(bName_).force() + J_deltaF * alpha_ / impact_dt_);
}

} // namespace mc_impact
