#include "ImpactAwareZMPConstraint.h"

namespace mc_impact
{

ImpactAwareZMPConstraint::ImpactAwareZMPConstraint(mi_qpEstimator & predictor,
                                                   const mc_rbdyn::Robot & robot,
                                                   const ImpactAwareConstraintParams<Eigen::Vector2d> & params)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), predictor_(predictor), robot_(robot),
  params_(params)
{

  if(enabledMcZMPArea())
  {
    mcZMPAreaPtr_ = std::make_shared<mc_impact::McZMPArea<Eigen::Vector2d>>(robot_, getParams().contactSetPtr,
                                                                            getParams().mcProjectionParams);
  }

  // Initialize the com Jacobian
  comJacobianPtr_ = std::make_shared<rbd::CoMJacobian>(robot_.mb());

  // Initialize the building blocks for the constraints.
  fixedSupportPolygonSetup_();

  std::cout << red << "Created ImpactAwareZMPConstraint." << reset << std::endl;
}

void ImpactAwareZMPConstraint::fixedSupportPolygonSetup_()
{
  // int numVertexZMP = mcZMPAreaPtr_->getNumVertex();
  int numVertexZMP = static_cast<int>(getParams().zmpAreaVertexSet.size());

  A_zmp_ = Eigen::MatrixXd::Zero(numVertexZMP, 6);

  // Modify the ieqConstraintBlocks.
  pointsToInequalityMatrix(getParams().zmpAreaVertexSet, ieqConstraintBlocksZMP_.G, ieqConstraintBlocksZMP_.h,
                           getParams().lowerSlope, getParams().upperSlope);
  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  A_zmp_.block(0, 0, numVertexZMP, 1) = getIeqBlocksZMP().G.block(0, 1, numVertexZMP, 1);
  /// A(:,1) = -G_x
  A_zmp_.block(0, 1, numVertexZMP, 1) = -getIeqBlocksZMP().G.block(0, 0, numVertexZMP, 1);
  /// A(:,5) = h
  A_zmp_.block(0, 5, numVertexZMP, 1) = -getIeqBlocksZMP().h;

  // Initialie the building blocks
  robotJointVelocity_.resize(dof_());

  std::cout << red << "Fixed-SupportPolygon set for the ImpactAwareZMPConstraint." << reset << std::endl;
}

void ImpactAwareZMPConstraint::updateMcAreas_(double height)
{
  if(getParams().updateMcZMPArea)
  {
    updateMcZMPAreas_(height);
    updateAbMatrixsize_();

    if(getParams().debug)
    {
      getMcZMPArea()->print();
    }
  }
}

void ImpactAwareZMPConstraint::updateMcZMPAreas_(double height)
{

  assert(mcZMPAreaPtr_ != nullptr);

  // --------------- (1) Update the Multi-contact ZMP area:

  mcZMPAreaPtr_->updateMcZMPArea(height);
  // int numVertex = static_cast<int>(iniVertexSet_.size());
  int numVertexZMP = getMcZMPArea()->getNumVertex();

  A_zmp_ = Eigen::MatrixXd::Zero(numVertexZMP, 6);

  // Set the inequality matrix blocks
  setIeqBlocksZMP(getMcZMPArea()->getIeqConstraint());

  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  A_zmp_.block(0, 0, numVertexZMP, 1) = getIeqBlocksZMP().G.block(0, 1, numVertexZMP, 1);
  /// A(:,1) = -G_x
  A_zmp_.block(0, 1, numVertexZMP, 1) = -getIeqBlocksZMP().G.block(0, 0, numVertexZMP, 1);
  /// A(:,5) = h
  A_zmp_.block(0, 5, numVertexZMP, 1) = -getIeqBlocksZMP().h;
}

void ImpactAwareZMPConstraint::getZMPBlocks(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench)
{
  exWrench.setZero();
  // int dof = predictor_.getSimRobot().mb().nrDof();
  // int dof = robot().mb().nrDof();
  sumJac.resize(6, dof_());
  sumJac.setZero();
  // sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  // (1) Go through the bodies with contact

  // for(auto idx = getParams().contactSetPtr->getContactMap().begin(); idx != getParams().contacts.end(); ++idx)

  for(auto & contactPair : getParams().contactSetPtr->getContactMap())
  {
    std::string bodyName = contactPair.second.getContactParams().bodyName;
    // sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(contactPair.second.getContactParams().bodyName).inv();
    sva::PTransformd X_ee_0 = robot().bodyPosW(contactPair.second.getContactParams().bodyName).inv();
    // sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(idx->bodyName).inv();
    // exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().forceSensor(idx->sensorName).wrench().vector();
    // exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(bodyName).vector();
    exWrench += X_ee_0.dualMatrix() * robot().bodyWrench(bodyName).vector();

    // sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(idx->bodyName);

    sumJac.block(0, 0, 3, dof_()) += X_ee_0.dualMatrix().block(0, 3, 3, 3) * predictor_.getJacobianDeltaF(bodyName);
    sumJac.block(3, 0, 3, dof_()) += X_ee_0.dualMatrix().block(3, 3, 3, 3) * predictor_.getJacobianDeltaF(bodyName);
    //
  } // end of for

  // (2) Go through the impacts
  if(getParams().multiContactCase)
  {
    for(auto impactIdx = predictor_.getImpactModels().begin(); impactIdx != predictor_.getImpactModels().end();
        ++impactIdx)
    {
      // sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();
      sva::PTransformd X_ee_0 = robot().bodyPosW(impactIdx->second->getImpactBody()).inv();

      sumJac.block(0, 0, 3, dof_()) +=
          X_ee_0.dualMatrix().block(0, 3, 3, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());

      sumJac.block(3, 0, 3, dof_()) +=
          X_ee_0.dualMatrix().block(3, 3, 3, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());

      // Add the hand force sensor measurement
      // exWrench += X_ee_0.dualMatrix() *
      // predictor_.getSimRobot().bodyWrench(impactIdx->second->getImpactBody()).vector();
    }
  }
}

void ImpactAwareZMPConstraint::calculateZMP_()
{

  // We use the joint velocity used by the impact model
  Eigen::VectorXd jointVel = predictor_.getImpactModels().begin()->second->getJointVel();
  double inv_t = (1 / predictor_.getImpactModels().begin()->second->getImpactDuration());

  // Eigen::Vector6d local_exWrench;
  // local_exWrench.setZero();

  Eigen::Vector6d contactBodyWrench;
  Eigen::Vector6d impactBodyWrench;
  contactBodyWrench.setZero();
  impactBodyWrench.setZero();

  Eigen::MatrixXd contactBodyWrenchJac; // Forces due to established contacts.
  Eigen::MatrixXd impactBodyWrenchJac; // Predicted forces due to expected impacts.

  contactBodyWrenchJac.resize(6, dof_());
  contactBodyWrenchJac.setZero();

  impactBodyWrenchJac.resize(6, dof_());
  impactBodyWrenchJac.setZero();

  // (1) Go through all the contact
  // Collect all the local Jacobians from all the contacts
  for(auto & contactPair : getParams().contactSetPtr->getContactMap())
  {
    std::string bodyName = contactPair.second.getContactParams().bodyName;
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(bodyName).inv();

    // Use condition "inContact" to exclude, e.g. hand contact.
    if(contactPair.second.inContact())
    {
      // contact established
      contactBodyWrench += X_ee_0.dualMatrix() * robot().bodyWrench(bodyName).vector();
      contactBodyWrenchJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(bodyName);
    }
    else
    {
      // contact not established, impact is expected.
      impactBodyWrench += X_ee_0.dualMatrix() * robot().bodyWrench(bodyName).vector();
      impactBodyWrenchJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(bodyName);
    }
  } // end of for

  // Sum of the impulse forces at the established contacts:
  Eigen::Vector6d impulseContactBodyWrench;
  impulseContactBodyWrench.setZero();
  impulseContactBodyWrench = inv_t * contactBodyWrenchJac * jointVel;
  Eigen::Vector6d contactBodyWrenchSum = impulseContactBodyWrench + contactBodyWrench;

  sva::ForceVecd sva_contactBodyWrenchSum;
  sva_contactBodyWrenchSum.force() = contactBodyWrenchSum.segment(3, 3); // Force starts at 3.
  sva_contactBodyWrenchSum.couple() = contactBodyWrenchSum.segment(0, 3); // Torque starts at 0.

  floatingBaseStates_.bipedalZMP.oneStepPreview =
      McZMPArea<Eigen::Vector2d>::zmpCalculation(Eigen::Vector3d::UnitZ(), sva_contactBodyWrenchSum);

  if(getParams().updateMcZMPArea)
  {
    floatingBaseStates_.bipedalZMP.current = getMcZMPArea()->getBipedalZMP();
  }
  else
  {
    sva::ForceVecd sva_ContactBodyWrench;
    sva_ContactBodyWrench.force() = contactBodyWrench.segment(3, 3); // Force starts at 3.
    sva_ContactBodyWrench.couple() = contactBodyWrench.segment(0, 3); // Torque starts at 0.

    floatingBaseStates_.bipedalZMP.current =
        McZMPArea<Eigen::Vector2d>::zmpCalculation(Eigen::Vector3d::UnitZ(), sva_ContactBodyWrench);
  }

  floatingBaseStates_.bipedalZMP.stateJump =
      floatingBaseStates_.bipedalZMP.oneStepPreview - floatingBaseStates_.bipedalZMP.current;

  Eigen::Vector6d sumWrench;
  sumWrench.setZero();
  // Add the predicted impulse
  sumWrench = inv_t * impactBodyWrenchJac * jointVel + contactBodyWrenchSum;

  sva::ForceVecd sva_sum;
  sva_sum.force() = sumWrench.segment(3, 3); // Force starts at 3.
  sva_sum.couple() = sumWrench.segment(0, 3); // Torque starts at 0.

  // zmpCalculation is a static function
  floatingBaseStates_.mcZMP.oneStepPreview =
      McZMPArea<Eigen::Vector2d>::zmpCalculation(Eigen::Vector3d::UnitZ(), sva_sum);

  if(getParams().updateMcZMPArea)
  {
    floatingBaseStates_.mcZMP.current = getMcZMPArea()->getMcZMP();
  }
  else
  {

    sva::ForceVecd sva_allContactBodyWrench;
    sva_allContactBodyWrench.force() =
        impactBodyWrench.segment(3, 3) + contactBodyWrench.segment(3, 3); // Force starts at 3.
    sva_allContactBodyWrench.couple() =
        impactBodyWrench.segment(0, 3) + contactBodyWrench.segment(0, 3); // Torque starts at 0.
    floatingBaseStates_.mcZMP.current =
        McZMPArea<Eigen::Vector2d>::zmpCalculation(Eigen::Vector3d::UnitZ(), sva_allContactBodyWrench);
  }

  floatingBaseStates_.mcZMP.stateJump = floatingBaseStates_.mcZMP.oneStepPreview - floatingBaseStates_.mcZMP.current;
}

void ImpactAwareZMPConstraint::updateFloatingBaseState_()
{
  // Update omega for DCM calculation
  calcOmega_(robot().com().z());

  Eigen::MatrixXd comJacobian = comJacobianPtr_->jacobian(robot().mb(), robot().mbc());
  Eigen::MatrixXd dcmJacobian = (comJacobian * predictor_.getJacobianDeltaAlpha()).block(0, 0, 2, dof_());

  floatingBaseStates_.Com = robot().com();
  floatingBaseStates_.ComAcc = robot().comAcceleration();

  // (1) COM velocity

  floatingBaseStates_.ComVel.current = robot().comVelocity();
  floatingBaseStates_.ComVel.stateJump = comJacobian * predictor_.getJointVelJump();
  floatingBaseStates_.ComVel.oneStepPreview = floatingBaseStates_.ComVel.current + floatingBaseStates_.ComVel.stateJump;

  // (2) ZMP
  calculateZMP_();

  // (3) DCM
  floatingBaseStates_.DCM.current.segment(0, 2) =
      floatingBaseStates_.Com.segment(0, 2) + floatingBaseStates_.ComVel.current.segment(0, 2) / getOmega();
  floatingBaseStates_.DCM.stateJump.segment(0, 2) =
      comJacobian.block(0, 0, 2, dof_()) * predictor_.getJointVelJump() / getOmega();
  floatingBaseStates_.DCM.oneStepPreview.segment(0, 2) =
      floatingBaseStates_.DCM.current.segment(0, 2) + floatingBaseStates_.DCM.stateJump.segment(0, 2);
}

void ImpactAwareZMPConstraint::compute()
{

  robotJointVelocity_ = (rbd::dofToVector(robot().mb(), robot().mbc().alpha));

  // Update the Multi-contact areas
  updateMcAreas_(2.0);

  // Update the ZMP, ComVel and DCM states
  updateFloatingBaseState_();

  if(getParams().debug)
  {
    std::cout << red << "ImpactAwareZMPConstraint computed FloatingBase-States." << reset << std::endl;
  }

  updateZMPConstraint_();
}

void ImpactAwareZMPConstraint::updateZMPConstraint_()
{

  Eigen::MatrixXd sumJacZMP;
  Eigen::Vector6d sumWrenchZMP;

  getZMPBlocks(sumJacZMP, sumWrenchZMP);

  // We assume that we always start from ZMP constraint unless it is not used.
  A_.block(0, 0, nrInEq(), dof_()) = (getParams().dt / getParams().impactDuration) * A_zmp_ * sumJacZMP;

  // Read the robot joint velocities.
  // rbd::paramToVector(robot.mbc().alpha, robotJointVelocity_);
  b_.segment(0, nrInEq()) = -(A_zmp_ * sumWrenchZMP + A_zmp_ * sumJacZMP * robotJointVelocity_ / getParams().dt);
}

void ImpactAwareZMPConstraint::updateAbMatrixsize_()
{
  A_.resize(nrInEq(), dof_());
  b_.resize(nrInEq());
}
void ImpactAwareZMPConstraint::printZMPConstraintMatrix() const
{

  std::cerr << red << "ZMP A matrix is: " << std::endl
            << cyan << A_ << std::endl
            << red << std::endl
            << "b vector is: " << cyan << b_.transpose() << std::endl
            << red << "Intermediat G matrix is: " << cyan << getIeqBlocksZMP().G << std::endl
            << red << "Intermediat h vector is: " << cyan << getIeqBlocksZMP().h.transpose() << reset << std::endl;
}

} // namespace mc_impact
