#include "zmpWithImpulse.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

zmpWithImpulse::zmpWithImpulse(mi_impactPredictor & predictor,
                               const std::vector<mc_impact::supportContact> & supports,
                               double dt,
                               double impact_dt,
                               const ZMPArea & area,
			       bool debug 
			       )
: InequalityConstraint(predictor.getRobot().robotIndex()), predictor_(predictor), dt_(dt), impact_dt_(impact_dt), supports_(supports), debug_(debug)
{

  // bName_ = bodyName;
  // sName_ = sensorName;
  A_zmp_ = Eigen::MatrixXd::Zero(4, 6);
  A_zmp_f_ = Eigen::MatrixXd::Zero(4, 6);

  // - cy - max_x * fz
  A_zmp_(0, 1) = -1;
  A_zmp_(0, 5) = -area.max_x;

  // cy + min_x * fz
  A_zmp_(1, 1) = 1;
  A_zmp_(1, 5) = area.min_x;

  // cx - max_y * fz
  A_zmp_(2, 0) = 1;
  A_zmp_(2, 5) = -area.max_y;

  // - cx + min_y * fz
  A_zmp_(3, 0) = -1;
  A_zmp_(3, 5) = area.min_y;


  A_zmp_f_ = A_zmp_;
  A_zmp_f_(0, 1) = 1;
  A_zmp_f_(1, 1) = -1;
  A_zmp_f_(2, 0) = -1;
  A_zmp_f_(3, 0) = 1;

  int nDof = predictor_.getRobot().mb().nrDof();

  alpha_.resize(nDof);
  // This needs double check
  A_.resize(4, nDof);
  b_.resize(4);
}

void zmpWithImpulse::getComItems(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench)
{
  exWrench.setZero();
  int dof = predictor_.getRobot().mb().nrDof();
  sumJac.resize(6, dof);
  sumJac.setZero();
  //sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  // Go through the bodies with contact  
  
  for(auto idx = supports_.begin(); idx != supports_.end(); ++idx)
  {
    sva::PTransformd X_ee_0 = predictor_.getRobot().bodyPosW(idx->bodyName).inv();
    exWrench += X_ee_0.dualMatrix()*predictor_.getRobot().forceSensor(idx->sensorName).wrench().vector();
   sumJac.block(0, 0, 3, dof) += X_ee_0.dualMatrix().block(0, 3, 3, 3)*predictor_.getJacobianDeltaF(idx->bodyName);
    sumJac.block(3, 0, 3, dof) += X_ee_0.dualMatrix().block(3, 3, 3, 3)*predictor_.getJacobianDeltaF(idx->bodyName);
  }// end of for

  // Add the impact body force
  sva::PTransformd X_ee_0 = predictor_.getRobot().bodyPosW(predictor_.getImpactBody()).inv();
  sumJac.block(0, 0, 3, dof) += X_ee_0.dualMatrix().block(0, 3, 3, 3)*predictor_.getJacobianDeltaF(predictor_.getImpactBody());
  sumJac.block(3, 0, 3, dof) += X_ee_0.dualMatrix().block(3, 3, 3, 3)*predictor_.getJacobianDeltaF(predictor_.getImpactBody());

}

void zmpWithImpulse::calcZMP_(const Eigen::MatrixXd & sumJac, const Eigen::Vector6d & exWrench)
{
 Eigen::VectorXd temp =
        (rbd::dofToVector(predictor_.getRobot().mb(), predictor_.getRobot().mbc().alpha)
         + rbd::dofToVector(predictor_.getRobot().mb(), predictor_.getRobot().mbc().alphaD) * predictor_.getImpactDuration_());

 Eigen::Vector6d impulseForceSum = (1/predictor_.getImpactDuration_())*sumJac*temp;
 double denominator = exWrench(5)  + impulseForceSum(5);

 zmpSensor_.x() = - (exWrench(1))/denominator;
 zmpSensor_.y() = (exWrench(0))/denominator;

 zmpPerturbation_.x() = - (impulseForceSum(1))/denominator;
 zmpPerturbation_.y() = (impulseForceSum(0))/denominator;
 
 zmpPrediction_ = zmpSensor_ + zmpPerturbation_;

}
void zmpWithImpulse::computeAb()
{
  const auto & robot = predictor_.getRobot();
  Eigen::MatrixXd sumJac;
  Eigen::Vector6d sumWrench;
  getComItems(sumJac, sumWrench);

  
  A_ = (dt_ / impact_dt_) * A_zmp_f_ * sumJac;
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  b_ = -(A_zmp_*sumWrench 
         + A_zmp_f_*sumJac* alpha_ / impact_dt_);
  if(debug_)
  {
    calcZMP_(sumJac, sumWrench);
    difference_ = A_* rbd::dofToVector(predictor_.getRobot().mb(), predictor_.getRobot().mbc().alphaD) - b_;
  }
}

} // namespace mc_impact
