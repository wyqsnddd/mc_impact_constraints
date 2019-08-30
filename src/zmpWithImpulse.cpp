/* Copyright 2018-2019 CNRS-UM LIRMM
 *
 * \author Yuquan Wang, Arnaud Tanguy and Pierre Gergondet
 *
 * 
 *
 * mc_impact_constraints is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * mc_impact_constraints is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with mc_impact_constraints. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "zmpWithImpulse.h"

namespace mc_impact
{

zmpWithImpulse::zmpWithImpulse(mi_qpEstimator & predictor,
                               const std::vector<mc_impact::supportContact> & supports,
                               double dt,
                               double impact_dt,
                               const ZMPArea & area,
                               bool allforce,
                               bool debug)
: InequalityConstraint(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt), impact_dt_(impact_dt),
  supports_(supports), area_(area), allForce_(allforce), debug_(debug)
{

  // bName_ = bodyName;
  // sName_ = sensorName;
  A_zmp_ = Eigen::MatrixXd::Zero(4, 6);
  A_zmp_f_ = Eigen::MatrixXd::Zero(4, 6);

  // - cy - max_x * fz
  A_zmp_(0, 1) = -1;
  A_zmp_(0, 5) = -area_.max_x;

  // cy + min_x * fz
  A_zmp_(1, 1) = 1;
  A_zmp_(1, 5) = area_.min_x;

  // cx - max_y * fz
  A_zmp_(2, 0) = 1;
  A_zmp_(2, 5) = -area_.max_y;

  // - cx + min_y * fz
  A_zmp_(3, 0) = -1;
  A_zmp_(3, 5) = area_.min_y;

  A_zmp_f_ = A_zmp_;
  A_zmp_f_(0, 1) = 1;
  A_zmp_f_(1, 1) = -1;
  A_zmp_f_(2, 0) = -1;
  A_zmp_f_(3, 0) = 1;

  int nDof = predictor_.getSimRobot().mb().nrDof();

  alpha_.resize(nDof);
  // This needs double check
  A_.resize(4, nDof);
  b_.resize(4);
}

void zmpWithImpulse::getInertialItems(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench)
{
  exWrench.setZero();
  int dof = predictor_.getSimRobot().mb().nrDof();
  sumJac.resize(6, dof);
  sumJac.setZero();
  // sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  // (1) Go through the bodies with contact

  for(auto idx = supports_.begin(); idx != supports_.end(); ++idx)
  {
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(idx->bodyName).inv();
    // exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().forceSensor(idx->sensorName).wrench().vector();
    exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(idx->bodyName).vector();

    //sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(idx->bodyName);

    
    sumJac.block(0, 0, 3, dof) += X_ee_0.dualMatrix().block(0, 3, 3, 3) * predictor_.getJacobianDeltaF(idx->bodyName);
    sumJac.block(3, 0, 3, dof) += X_ee_0.dualMatrix().block(3, 3, 3, 3) * predictor_.getJacobianDeltaF(idx->bodyName);
    //
  } // end of for

  // (2) Go through the impacts
  if(allForce_)
  {
    for(auto impactIdx = predictor_.getImpactModels().begin(); impactIdx != predictor_.getImpactModels().end();
        ++impactIdx)
    {
      sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();

      /*
      sumJac +=
          X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());

	  */
      
       sumJac.block(0, 0, 3, dof) +=
           X_ee_0.dualMatrix().block(0, 3, 3, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());

       sumJac.block(3, 0, 3, dof) +=
           X_ee_0.dualMatrix().block(3, 3, 3, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());
     

       // Add the hand force sensor measurement 
       //exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(impactIdx->second->getImpactBody()).vector();

    }
  }
}

void zmpWithImpulse::calcZMP_()
{
/*
  Eigen::VectorXd temp = (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alpha)
                          + rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD)
                                * predictor_.getImpactModels().begin()->second->getTimeStep());
				* */
  Eigen::VectorXd temp = predictor_.getImpactModels().begin()->second->getJointVel();
  double inv_t = (1 / predictor_.getImpactModels().begin()->second->getImpactDuration());

  Eigen::Vector6d local_exWrench;
  local_exWrench.setZero();
  int dof = predictor_.getSimRobot().mb().nrDof();
  Eigen::MatrixXd local_sumJac;
  local_sumJac.resize(6, dof);
  local_sumJac.setZero();
  // sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  // (1) Go through the bodies with contact

  for(auto idx = supports_.begin(); idx != supports_.end(); ++idx)
  {
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(idx->bodyName).inv();
    local_exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(idx->bodyName).vector();
    local_sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(idx->bodyName);

  } // end of for

  Eigen::Vector6d impulseForceSum;
  impulseForceSum.setZero();
  impulseForceSum = inv_t * local_sumJac * temp;
  double denominator = local_exWrench(5) + impulseForceSum(5);

  zmpPrediction_feetforce_.x() = -(local_exWrench(1) + impulseForceSum(1)) / denominator;
  zmpPrediction_feetforce_.y() = (local_exWrench(0) + impulseForceSum(0)) / denominator;

  // (2) Go through the impacts
  for(auto impactIdx = predictor_.getImpactModels().begin(); impactIdx != predictor_.getImpactModels().end();
      ++impactIdx)
  {
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();

    local_sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());
  }
  // Recalculate impulseSum:
  impulseForceSum = inv_t * local_sumJac * temp;
  denominator = local_exWrench(5) + impulseForceSum(5);

  zmpSensor_.x() = -(local_exWrench(1)) / denominator;
  zmpSensor_.y() = (local_exWrench(0)) / denominator;

  zmpPerturbation_.x() = -(impulseForceSum(1)) / denominator;
  zmpPerturbation_.y() = (impulseForceSum(0)) / denominator;

  zmpPrediction_allforce_ = zmpSensor_ + zmpPerturbation_;
}

void zmpWithImpulse::computeAb()
{
  const auto & robot = predictor_.getSimRobot();
  Eigen::MatrixXd sumJac;
  Eigen::Vector6d sumWrench;
  getInertialItems(sumJac, sumWrench);

  //A_ = (dt_ / impact_dt_) * A_zmp_f_ * sumJac;
  A_ = (dt_ / impact_dt_) * A_zmp_ * sumJac;
  /*
    alpha_ =
          (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alpha));
    */
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  //b_ = -(A_zmp_ * sumWrench + A_zmp_f_ * sumJac * alpha_ / impact_dt_);
  b_ = -(A_zmp_ * sumWrench + A_zmp_ * sumJac * alpha_ / impact_dt_);

  if(debug_)
  {
    calcZMP_();
    difference_ = A_ * rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD) - b_;
  }
}

} // namespace mc_impact
