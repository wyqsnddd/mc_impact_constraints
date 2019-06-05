#include "zmpWithImpulse.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

zmpWithImpulse::zmpWithImpulse(mi_impactPredictor & predictor,
                               // const std::string & bodyName,
                               // const std::string & sensorName,
                               const std::vector<mc_impact::supportContact> & supports,
                               double dt,
                               double impact_dt,
                               const ZMPArea & area)
: InequalityConstraint(predictor.getRobot().robotIndex()), predictor_(predictor), dt_(dt), impact_dt_(impact_dt),
  supports_(supports)
{

  // bName_ = bodyName;
  // sName_ = sensorName;

  A_zmp_ = Eigen::MatrixXd::Zero(4, 6);

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
  sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  
  for(auto idx = supports_.begin(); idx != supports_.end(); ++idx)
  {
    sva::PTransformd X_0_ee = predictor_.getRobot().bodyPosW(idx->bodyName);
    sva::PTransformd X_ee_CoM = X_0_CoM * X_0_ee.inv();
    //Eigen::Matrix6d tMatrix = X_ee_CoM.dualMatrix();
    exWrench += X_ee_CoM.dualMatrix()*predictor_.getRobot().forceSensor(idx->sensorName).wrench().vector();
    sumJac.block(0, 0, 3, dof) += X_ee_CoM.dualMatrix().block(0, 3, 3, 3)*predictor_.getJacobianDeltaF(idx->bodyName);
    sumJac.block(3, 0, 3, dof) += X_ee_CoM.dualMatrix().block(3, 3, 3, 3)*predictor_.getJacobianDeltaF(idx->bodyName);
  }// end of for
}
void zmpWithImpulse::computeAb()
{

  const auto & robot = predictor_.getRobot();
  Eigen::MatrixXd sumJac;
  Eigen::Vector6d sumWrench;
  getComItems(sumJac, sumWrench);
 
  //const auto & J_deltaF = predictor_.getJacobianDeltaF(bName_);
  // std::cout<<"size of J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;
  // std::cout<<"size of reduced J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;

  // A_ = (dt_ / impact_dt_) * A_zmp_.block(0, 3, 4, 3)*J_deltaF.block(0, startIndex_, J_deltaF.rows(), J_deltaF.cols()
  // - startIndex_);
  A_ = (dt_ / impact_dt_) * A_zmp_ * sumJac;
  // std::cout<<"size of A_: "<<A_.rows()<<", "<<A_.cols()<<std::endl;
  rbd::paramToVector(robot.mbc().alpha, alpha_);

  // std::cout<<"size of alpha_"<<alpha_.rows()<<std::endl;
  b_ = -A_zmp_*(sumWrench 
         + sumJac* alpha_ / impact_dt_);
}

} // namespace mc_impact
