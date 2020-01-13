#include "ZMPWithImpulse.h"

namespace mc_impact
{

template<typename Point>
ZMPWithImpulse<Point>::ZMPWithImpulse(mi_qpEstimator & predictor,

                                      std::shared_ptr<mc_impact::McZMPArea<Point>> mcZMPAreaPtr,
                                      const ImpactAwareConstraintParams<Point> & params)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), predictor_(predictor),
  mcZMPAreaPtr_(mcZMPAreaPtr), params_(params)
{
  int numVertex = static_cast<int>(getParams().zmpAreaVertexSet.size());
  A_zmp_ = Eigen::MatrixXd::Zero(numVertex, 6);

  // Modify the ieqConstraintBlocks. 
  pointsToInequalityMatrix(getParams().zmpAreaVertexSet, ieqConstraintBlocks_.G_zmp, ieqConstraintBlocks_.h_zmp, 
                                   getParams().lowerSlope, getParams().upperSlope);


  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  A_zmp_.block(0, 0, numVertex, 1) = getIeqBlocks().G_zmp.block(0, 1, numVertex, 1);
  /// A(:,1) = -G_x
  A_zmp_.block(0, 1, numVertex, 1) = -getIeqBlocks().G_zmp.block(0, 0, numVertex, 1);
  /// A(:,5) = h
  A_zmp_.block(0, 5, numVertex, 1) = -getIeqBlocks().h_zmp;

  int nDof = predictor_.getSimRobot().mb().nrDof();

  alpha_.resize(nDof);
  A_.resize(numVertex, nDof);
  b_.resize(numVertex);

}

template<typename Point>
void ZMPWithImpulse<Point>::computeMcZMPArea_(double height)
{

  // Update the Multi-contact ZMP area. 
  getMcZMPArea()->computeMcZMPArea(height);

  // int numVertex = static_cast<int>(iniVertexSet_.size());
  int numVertex = getMcZMPArea()->getNumVertex();

  /*
  A_zmp_ = Eigen::MatrixXd::Zero(numVertex, 6);
  // Set the inequality matrix blocks
  setIeqBlocks(getMcZMPArea()->getIeqConstraint());

  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  A_zmp_.block(0, 0, numVertex, 1) = getIeqBlocks().G_zmp.block(0, 1, numVertex, 1);
  /// A(:,1) = -G_x
  A_zmp_.block(0, 1, numVertex, 1) = -getIeqBlocks().G_zmp.block(0, 0, numVertex, 1);
  /// A(:,5) = h
  A_zmp_.block(0, 5, numVertex, 1) = -getIeqBlocks().h_zmp;

  */
}

template<typename Point>
void ZMPWithImpulse<Point>::getInertialItems(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench)
{
  exWrench.setZero();
  int dof = predictor_.getSimRobot().mb().nrDof();
  sumJac.resize(6, dof);
  sumJac.setZero();
  // sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  // (1) Go through the bodies with contact

  // for(auto idx = getParams().contactSetPtr->getContactMap().begin(); idx != getParams().contacts.end(); ++idx)

  for(auto & contactPair : getParams().contactSetPtr->getContactMap())
  {
    std::string bodyName = contactPair.second.getContactParams().bodyName;
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(contactPair.second.getContactParams().bodyName).inv();
    // sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(idx->bodyName).inv();
    // exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().forceSensor(idx->sensorName).wrench().vector();
    exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(bodyName).vector();

    // sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(idx->bodyName);

    sumJac.block(0, 0, 3, dof) += X_ee_0.dualMatrix().block(0, 3, 3, 3) * predictor_.getJacobianDeltaF(bodyName);
    sumJac.block(3, 0, 3, dof) += X_ee_0.dualMatrix().block(3, 3, 3, 3) * predictor_.getJacobianDeltaF(bodyName);
    //
  } // end of for

  // (2) Go through the impacts
  if(getParams().multiContactCase)
  {
    for(auto impactIdx = predictor_.getImpactModels().begin(); impactIdx != predictor_.getImpactModels().end();
        ++impactIdx)
    {
      sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();

      sumJac.block(0, 0, 3, dof) +=
          X_ee_0.dualMatrix().block(0, 3, 3, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());

      sumJac.block(3, 0, 3, dof) +=
          X_ee_0.dualMatrix().block(3, 3, 3, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());

      // Add the hand force sensor measurement
      // exWrench += X_ee_0.dualMatrix() *
      // predictor_.getSimRobot().bodyWrench(impactIdx->second->getImpactBody()).vector();
    }
  }
}

template<typename Point>
void ZMPWithImpulse<Point>::calcZMP_()
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

  // for(auto idx = getParams().contactSetPtr->getContactMap().begin(); idx != getParams().contacts.end(); ++idx)
  for(auto & contactPair : getParams().contactSetPtr->getContactMap())
  {
    std::string bodyName = contactPair.second.getContactParams().bodyName;
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(bodyName).inv();
    local_exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(bodyName).vector();
    local_sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(bodyName);

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

    local_sumJac +=
        X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());
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

template<typename Point>
bool ZMPWithImpulse<Point>::pointInsideSupportPolygon(const Point & input)
{

  Eigen::VectorXd result = getIeqBlocks().G_zmp * input - getIeqBlocks().h_zmp;
  if(getParams().debug)
  {
    std::cerr<<"G_zmp: is: "<<std::endl<<getIeqBlocks().G_zmp<<std::endl;
    std::cerr<<"h_zmp: is: "<<std::endl<<getIeqBlocks().h_zmp.transpose()<<std::endl;

    std::cerr<<"The test result is: "<< result.transpose()<<std::endl; 
    std::cerr<<"The zmp area vertices size is: "<< getParams().zmpAreaVertexSet.size()<<std::endl;
    
  }
for(int ii = 0; ii < static_cast<int>(getParams().zmpAreaVertexSet.size()); ii++)
    {
      if(result(ii) > 0) return false;
    }
  return true;
}

template<typename Point>
void ZMPWithImpulse<Point>::compute()
{
  const auto & robot = predictor_.getSimRobot();
  Eigen::MatrixXd sumJac;
  Eigen::Vector6d sumWrench;
  getInertialItems(sumJac, sumWrench);

  if(updateMcZMPArea())
  {
    computeMcZMPArea_(2.0);
  }

  A_ = (getParams().dt / getParams().impactDuration) * A_zmp_ * sumJac;
  /*
    alpha_ =
          (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alpha));
  */

  rbd::paramToVector(robot.mbc().alpha, alpha_);
  b_ = -(A_zmp_ * sumWrench + A_zmp_ * sumJac * alpha_ / getParams().dt);

  if(getParams().debug)
  {
    calcZMP_();
    difference_ = A_ * rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD) - b_;
  }
}

// The explicit instantiation
template struct mc_impact::ZMPWithImpulse<Eigen::Vector2d>;

} // namespace mc_impact
