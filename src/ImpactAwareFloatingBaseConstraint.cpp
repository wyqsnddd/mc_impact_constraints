#include "ImpactAwareFloatingBaseConstraint.h"

namespace mc_impact
{
ImpactAwareFloatingBaseConstraint::ImpactAwareFloatingBaseConstraint(
		std::shared_ptr<mi_qpEstimator> predictorPtr,
    std::shared_ptr<McContactSet> contactSetPtr,
    const ImpactAwareConstraintParams<Eigen::Vector2d> & params
    )
: mc_solver::InequalityConstraintRobot(getPredictor()->getSimRobot().robotIndex()), predictorPtr_(predictorPtr), params_(params), contactSetPtr_(contactSetPtr), robot_(getPredictor()->getSimRobot())
{

  // Set the constraining status of the floating-base state.
  constrainingStatus_ = 0;
  if(getParams().constrainingZMP)
  {
    constrainingStatus_ += 1;
  }
  if(getParams().constrainingDCM)
  {
    constrainingStatus_ += 3;
  }

  // When the constraints are enabled and  the realtime computation of McComArea and McDCMArea are needed.
  // As long as McZMPArea or McDCMArea need to be udpated
  // if((getParams().updateMcZMPArea or getParams().updateMcDCMArea))
  if(enabledMcZMPArea())
  {

    mcZMPAreaPtr_ = std::make_shared<mc_impact::McZMPArea<Eigen::Vector2d>>(robot_, contactSetPtr_,
                                                                            getParams().mcProjectionParams);
    // mcZMPAreaPtr_ = std::make_shared<mc_impact::McZMPArea<Eigen::Vector2d>>(getPredictor()->.getSimRobot(),
    // getParams().contactSetPtr, getParams().mcProjectionParams);
  }

  // When the constraints are enabled and  the realtime computation of McComArea and McDCMArea are needed.
  // if(getParams().constrainingDCM and getParams().updateMcDCMArea)
  if(enabledMcDCMArea())
  {

    assert(mcZMPAreaPtr_ != nullptr);

    /*
    mcComAreaPtr_ =
      std::make_shared<mc_impact::McComArea>(getPredictor()->.getSimRobot(), getParams().contactSetPtr,
    getParams().mcProjectionParams); mcDCMAreaPtr_ = std::make_shared<mc_impact::McDCMArea>(mcZMPAreaPtr_,
    mcComAreaPtr_);

    */

    mcComAreaPtr_ =
        std::make_shared<mc_impact::McComArea>(robot_, contactSetPtr_, getParams().mcProjectionParams);

    mcDCMAreaPtr_ = std::make_shared<mc_impact::McDCMArea>(mcZMPAreaPtr_, mcComAreaPtr_);
  }
  // Initialize the com Jacobian
  comJacobianPtr_ = std::make_shared<rbd::CoMJacobian>(robot_.mb());

  // Initialize the building blocks for the constraints.
  fixedSupportPolygonSetup_();

  std::cout << red << "Created ImpactAwareFloatingBaseConstraint." << reset << std::endl;
}

ImpactAwareFloatingBaseConstraint::~ImpactAwareFloatingBaseConstraint()
{

 if(checkLoggersAdded_())
 {
   removeLogFloatingBaseStates(); 
 }

}


void ImpactAwareFloatingBaseConstraint::fixedSupportPolygonSetup_()
{
  if(zmpConstraintEnabled())
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
  }
  // int numVertexDCM = mcDCMAreaPtr_->getNumVertex();
  // A_dcm_ = Eigen::MatrixXd::Zero(numVertexDCM, 6);
  // Modify the ieqConstraintBlocks.
  if(dcmConstraintEnabled())
  {
    pointsToInequalityMatrix(getParams().dcmAreaVertexSet, ieqConstraintBlocksDCM_.G, ieqConstraintBlocksDCM_.h,
                             getParams().lowerSlope, getParams().upperSlope);

    /// Needs to be checked carefully, compare to the 4 dim case
    // A(:,0) = G_y
    // A_dcm_.block(0, 0, numVertexDCM, 1) = getIeqBlocksDCM().G.block(0, 1, numVertexDCM, 1);
    /// A(:,1) = -G_x
    // A_dcm_.block(0, 1, numVertexDCM, 1) = -getIeqBlocksDCM().G.block(0, 0, numVertexDCM, 1);
    /// A(:,5) = h
    // A_dcm_.block(0, 5, numVertexDCM, 1) = -getIeqBlocksDCM().h;
  }

  // Initialie the building blocks
  robotJointVelocity_.resize(dof_());

  /*
  A_.resize(nrInEq(), dof_());
  b_.resize(nrInEq());
  */

  std::cout << red << "Fixed-SupportPolygon set for the ImpactAwareFloatingBaseConstraint." << reset << std::endl;
}
void ImpactAwareFloatingBaseConstraint::updateMcZMPAreas_(double height)
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

void ImpactAwareFloatingBaseConstraint::updateMcAreas_(double height)
{
  switch(getConstrainingStatus())
  {
    case 1:
      if(getParams().updateMcZMPArea)
      {
        updateMcZMPAreas_(height);
        if(getParams().debug)
        {
          getMcZMPArea()->print();
        }
      }
      break;
    case 3:
      if(getParams().updateMcDCMArea)
      {
        updateMcDCMAreas_();
        if(getParams().debug)
        {
          getMcZMPArea()->print();
          getMcComArea()->print();
          getMcDCMArea()->print();
        }
      }
      break;
    case 4:

      if(getParams().updateMcZMPArea)
      {
        updateMcZMPAreas_(height);

        if(getParams().debug)
        {
          getMcZMPArea()->print();
        }
      }
      if(getParams().updateMcDCMArea)
      {
        updateMcDCMAreas_();

        if(getParams().debug)
        {
          getMcZMPArea()->print();
          getMcComArea()->print();
          getMcDCMArea()->print();
        }
      }
      break;
    default:
      throw std::runtime_error("The constraining status is not set correctly");
  }
}
void ImpactAwareFloatingBaseConstraint::updateMcDCMAreas_()
{
  assert(mcZMPAreaPtr_ != nullptr);
  assert(mcComAreaPtr_ != nullptr);
  assert(mcDCMAreaPtr_ != nullptr);

  // --------------- (2) Update the Multi-contact ZMP area:

  // If ZMP constraint is not used, we also need to update the ZMP area:

  if(not zmpConstraintEnabled())
  {
    mcZMPAreaPtr_->updateMcZMPArea(2.0);
  }

  mcComAreaPtr_->updateMcComArea();
  mcDCMAreaPtr_->updateMcDCMArea();

  // Set the inequality matrix blocks
  setIeqBlocksDCM(getMcDCMArea()->getIeqConstraint());
}

void ImpactAwareFloatingBaseConstraint::getZMPBlocks(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench)
{
  exWrench.setZero();
  // int dof = getPredictor()->.getSimRobot().mb().nrDof();
  // int dof = robot().mb().nrDof();
  sumJac.resize(6, dof_());
  sumJac.setZero();
  // sva::PTransformd X_0_CoM = sva::PTransformd(getPredictor()->.getRobot().com());
  // (1) Go through the bodies with contact

  // for(auto idx = getParams().contactSetPtr->getContactMap().begin(); idx != getParams().contacts.end(); ++idx)

  for(auto & contactPair : contactSetPtr_->getContactMap())
  {
    std::string bodyName = contactPair.second.getContactParams().bodyName;
    // sva::PTransformd X_ee_0 = getPredictor()->.getSimRobot().bodyPosW(contactPair.second.getContactParams().bodyName).inv();
    sva::PTransformd X_ee_0 = robot().bodyPosW(contactPair.second.getContactParams().bodyName).inv();
    // sva::PTransformd X_ee_0 = getPredictor()->.getSimRobot().bodyPosW(idx->bodyName).inv();
    // exWrench += X_ee_0.dualMatrix() * getPredictor()->.getSimRobot().forceSensor(idx->sensorName).wrench().vector();
    // exWrench += X_ee_0.dualMatrix() * getPredictor()->.getSimRobot().bodyWrench(bodyName).vector();
    exWrench += X_ee_0.dualMatrix() * robot().bodyWrench(bodyName).vector();

    // sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * getPredictor()->.getJacobianDeltaF(idx->bodyName);

    sumJac.block(0, 0, 3, dof_()) += X_ee_0.dualMatrix().block(0, 3, 3, 3) * getPredictor()->getJacobianDeltaF(bodyName);
    sumJac.block(3, 0, 3, dof_()) += X_ee_0.dualMatrix().block(3, 3, 3, 3) * getPredictor()->getJacobianDeltaF(bodyName);
    //
  } // end of for

  // (2) Go through the impacts
  if(getParams().multiContactCase)
  {
    for(auto impactIdx = getPredictor()->getImpactModels().begin(); impactIdx != getPredictor()->getImpactModels().end();
        ++impactIdx)
    {
      // sva::PTransformd X_ee_0 = getPredictor()->.getSimRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();
      sva::PTransformd X_ee_0 = robot().bodyPosW(impactIdx->second->getImpactBody()).inv();

      sumJac.block(0, 0, 3, dof_()) +=
          X_ee_0.dualMatrix().block(0, 3, 3, 3) * getPredictor()->getJacobianDeltaF(impactIdx->second->getImpactBody());

      sumJac.block(3, 0, 3, dof_()) +=
          X_ee_0.dualMatrix().block(3, 3, 3, 3) * getPredictor()->getJacobianDeltaF(impactIdx->second->getImpactBody());

      // Add the hand force sensor measurement
      // exWrench += X_ee_0.dualMatrix() *
      // getPredictor()->.getSimRobot().bodyWrench(impactIdx->second->getImpactBody()).vector();
    }
  }
}

void ImpactAwareFloatingBaseConstraint::calculateZMP_()
{

  // We use the joint velocity used by the impact model
  Eigen::VectorXd jointVel = getPredictor()->getImpactModels().begin()->second->getJointVel();
  double inv_t = (1 / getPredictor()->getImpactModels().begin()->second->getImpactDuration());

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
  for(auto & contactPair : contactSetPtr_->getContactMap())
  {

    std::string bodyName = contactPair.second.getContactParams().bodyName;
    sva::PTransformd X_ee_0 = getPredictor()->getSimRobot().bodyPosW(bodyName).inv();

    // Use condition "inContact" to exclude, e.g. hand contact.
    if(contactPair.second.inContact())
    {
      // contact established
      contactBodyWrench += X_ee_0.dualMatrix() * robot().bodyWrench(bodyName).vector();
      contactBodyWrenchJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * getPredictor()->getJacobianDeltaF(bodyName);
    }
    else
    {
      // contact not established, impact is expected.
      impactBodyWrench += X_ee_0.dualMatrix() * robot().bodyWrench(bodyName).vector();
      impactBodyWrenchJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * getPredictor()->getJacobianDeltaF(bodyName);
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

  // zmpCalculation is a static function
  floatingBaseStates_.bipedalZMP.oneStepPreview =
      McZMPArea<Eigen::Vector2d>::zmpCalculation(Eigen::Vector3d::UnitZ(), sva_contactBodyWrenchSum);

  // TODO:  we need to calculate with the static function if realtime update is not enabled.
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

void ImpactAwareFloatingBaseConstraint::updateFloatingBaseState_()
{
  // Update omega for DCM calculation
  calcOmega(robot().com().z());

  Eigen::MatrixXd comJacobian = comJacobianPtr_->jacobian(robot().mb(), robot().mbc());
  Eigen::MatrixXd dcmJacobian = (comJacobian * getPredictor()->getJacobianDeltaAlpha()).block(0, 0, 2, dof_());

  floatingBaseStates_.Com = robot().com();
  floatingBaseStates_.ComAcc = robot().comAcceleration();

  // (1) COM velocity

  floatingBaseStates_.ComVel.current = robot().comVelocity();
  floatingBaseStates_.ComVel.stateJump = comJacobian * getPredictor()->getJointVelJump();
  floatingBaseStates_.ComVel.oneStepPreview = floatingBaseStates_.ComVel.current + floatingBaseStates_.ComVel.stateJump;

  // (2) ZMP
  calculateZMP_();

  // (3) DCM
  floatingBaseStates_.DCM.current.segment(0, 2) =
      floatingBaseStates_.Com.segment(0, 2) + floatingBaseStates_.ComVel.current.segment(0, 2) / getOmega();
  floatingBaseStates_.DCM.stateJump.segment(0, 2) =
      comJacobian.block(0, 0, 2, dof_()) * getPredictor()->getJointVelJump() / getOmega();
  floatingBaseStates_.DCM.oneStepPreview.segment(0, 2) =
      floatingBaseStates_.DCM.current.segment(0, 2) + floatingBaseStates_.DCM.stateJump.segment(0, 2);
}

void ImpactAwareFloatingBaseConstraint::updateDCMConstraint_()
{

  assert(getParams().constrainingDCM == true);
  // Calculates the new DCM constraint blocks internally.

  Eigen::MatrixXd jacDCM = (comJacobianPtr_->jacobian(robot().mb(), robot().mbc()) * getPredictor()->getJacobianDeltaAlpha())
                               .block(0, 0, 2, dof_());

  A_.block(getDCMRowNr_(), 0, getDCMConstraintSize_(), dof_()) =
      getParams().dt / getOmega() * ieqConstraintBlocksDCM_.G * jacDCM;

  b_.segment(getDCMRowNr_(), getDCMConstraintSize_()) =
      (ieqConstraintBlocksDCM_.h - ieqConstraintBlocksDCM_.G * floatingBaseStates_.DCM.current.segment(0, 2))
      - ieqConstraintBlocksDCM_.G * jacDCM * robotJointVelocity_ / getOmega();
}

void ImpactAwareFloatingBaseConstraint::updateZMPConstraint_()
{

  assert(getParams().constrainingZMP == true);

  Eigen::MatrixXd sumJacZMP;
  Eigen::Vector6d sumWrenchZMP;
  getZMPBlocks(sumJacZMP, sumWrenchZMP);

  // We assume that we always start from ZMP constraint unless it is not used.
  A_.block(0, 0, getZMPConstraintSize_(), dof_()) = (getParams().dt / getParams().impactDuration) * A_zmp_ * sumJacZMP;

  // Read the robot joint velocities.
  // rbd::paramToVector(robot.mbc().alpha, robotJointVelocity_);
  b_.segment(0, getZMPConstraintSize_()) =
      -(A_zmp_ * sumWrenchZMP + A_zmp_ * sumJacZMP * robotJointVelocity_ / getParams().dt);
}

void ImpactAwareFloatingBaseConstraint::compute()
{

  // std::cout<< cyan <<" This is a test sentence"<< reset<<std::endl;
  // const auto & robot = getPredictor()->.getSimRobot();
  // Read the SIMULATED robot joint velocity:

  robotJointVelocity_ = (rbd::dofToVector(robot().mb(), robot().mbc().alpha));

  // Update the Multi-contact areas
  updateMcAreas_(2.0);

  // Update the ZMP, ComVel and DCM states
  updateFloatingBaseState_();

  if(getParams().debug)
  {
    std::cout << red << "ImpactAwareFloatingBaseConstraint computed FloatingBase-States." << reset << std::endl;
  }
  // std::cout<<red <<"Updated McAreas"<< reset<<std::endl;

  // Update the Constraint size based on the number of the Multi-contact areas vertices.
  updateAbMatrixsize_();

  if(getParams().debug)
  {

    std::cout << red << "Updated McAreas" << reset << std::endl;
    std::cout << red << "ImpactAwareFloatingBaseConstraint status: " << getConstrainingStatus() << reset << std::endl;

    std::cout << red << "A_ size: " << A_.rows() << ", " << A_.cols() << reset << std::endl;
    std::cout << red << "b_ size: " << b_.rows() << std::endl;

    if(zmpConstraintEnabled())
    {
      printZMPConstraintMatrix();
    }
    if(dcmConstraintEnabled())
    {
      printDCMConstraintMatrix();
    }
  }

  // Update the Constraints formulation
  switch(getConstrainingStatus())
  {
    case 1:
      updateZMPConstraint_();
      break;
    case 3:
      updateDCMConstraint_();
      break;
    case 4:
      updateZMPConstraint_();
      updateDCMConstraint_();
      break;
    default:
      throw std::runtime_error("The constraining status is not set correctly");
  }
}

bool ImpactAwareFloatingBaseConstraint::pointInsideMcZMPArea(const Eigen::Vector3d & samplePoint) const
{

  Eigen::VectorXd result = getIeqBlocksZMP().G * samplePoint - getIeqBlocksZMP().h;

  for(int ii = 0; ii < static_cast<int>(result.size()); ii++)
  {
    if(result(ii) > 0) return false;
  }
  return true;
}

bool ImpactAwareFloatingBaseConstraint::pointInsideMcDCMArea(const Eigen::Vector3d & samplePoint) const
{

  Eigen::VectorXd result = getIeqBlocksDCM().G * samplePoint - getIeqBlocksDCM().h;

  for(int ii = 0; ii < static_cast<int>(result.size()); ii++)
  {
    if(result(ii) > 0) return false;
  }
  return true;
}

void ImpactAwareFloatingBaseConstraint::updateAbMatrixsize_()
{
  A_.resize(nrInEq(), dof_());
  b_.resize(nrInEq());
}

void ImpactAwareFloatingBaseConstraint::printZMPConstraintMatrix() const
{

  std::cerr << red << "ZMP A matrix is: " << std::endl
            << cyan << A_ << std::endl
            << red << std::endl
            << "b vector is: " << cyan << b_.transpose() << std::endl
            << red << "Intermediat G matrix is: " << cyan << getIeqBlocksZMP().G << std::endl
            << red << "Intermediat h vector is: " << cyan << getIeqBlocksZMP().h.transpose() << reset << std::endl;
}

void ImpactAwareFloatingBaseConstraint::printDCMConstraintMatrix() const
{

  std::cerr << red << "DCM A matrix is: " << std::endl
            << cyan << A_ << std::endl
            << red << std::endl
            << "b vector is: " << cyan << b_.transpose() << std::endl
            << red << "Intermediat G matrix is: " << cyan << getIeqBlocksDCM().G << std::endl
            << red << "Intermediat h vector is: " << cyan << b_.transpose() << reset << std::endl;
}

void ImpactAwareFloatingBaseConstraint::logFloatingBaseStates(mc_control::fsm::Controller & ctl) 
{

  setControllerWithLogger_(ctl);
  
  // DCM states
  // auto & ctl = static_cast<mc_impact::Controller &>(ctlInput);

  ctl.logger().addLogEntry("FloatingBaseState_DCM_Prediction",
                           [this]() { return this->getFloatingBaseStates().DCM.oneStepPreview; });

  ctl.logger().addLogEntry("FloatingBaseState_DCM_StateJump",
                           [this]() { return this->getFloatingBaseStates().DCM.stateJump; });

  ctl.logger().addLogEntry("FloatingBaseState_DCM", [this]() { return this->getFloatingBaseStates().DCM.current; });

  // Multi-contact ZMP states

  ctl.logger().addLogEntry("FloatingBaseState_ZMP_Mc_Prediction",
                           [this]() { return getFloatingBaseStates().mcZMP.oneStepPreview; });
  ctl.logger().addLogEntry("FloatingBaseState_ZMP_Mc_StateJump",
                           [this]() { return getFloatingBaseStates().mcZMP.stateJump; });
  ctl.logger().addLogEntry("FloatingBaseState_ZMP_Mc", [this]() { return getFloatingBaseStates().mcZMP.current; });

  // Bipedal ZMP states
  ctl.logger().addLogEntry("FloatingBaseState_ZMP_Bipedal_Prediction",
                           [this]() { return getFloatingBaseStates().bipedalZMP.oneStepPreview; });
  ctl.logger().addLogEntry("FloatingBaseState_ZMP_Bipedal_StateJump",
                           [this]() { return getFloatingBaseStates().bipedalZMP.stateJump; });
  ctl.logger().addLogEntry("FloatingBaseState_ZMP_Bipedal",
                           [this]() { return getFloatingBaseStates().bipedalZMP.current; });

  // Com states
  ctl.logger().addLogEntry("FloatingBaseState_Com", [this]() { return getFloatingBaseStates().Com; });
  ctl.logger().addLogEntry("FloatingBaseState_Com_Acc", [this]() { return getFloatingBaseStates().ComAcc; });

  // Com-vel states
  ctl.logger().addLogEntry("FloatingBaseState_Com_Vel_Prediction",
                           [this]() { return getFloatingBaseStates().ComVel.oneStepPreview; });
  ctl.logger().addLogEntry("FloatingBaseState_Com_Vel_StateJump",
                           [this]() { return getFloatingBaseStates().ComVel.stateJump; });
  ctl.logger().addLogEntry("FloatingBaseState_Com_Vel", [this]() { return getFloatingBaseStates().ComVel.current; });
}

void ImpactAwareFloatingBaseConstraint::removeLogFloatingBaseStates() 
{

  if(checkLoggersAdded_())
  {
    // DCM states
    // auto & ctl = static_cast<mc_impact::Controller &>(ctlInput);

    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_DCM_Prediction");

    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_DCM_StateJump");

    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_DCM");

    // Multi-contact ZMP states

    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_ZMP_Mc_Prediction");
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_ZMP_Mc");

    // Bipedal ZMP states
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_ZMP_Bipedal_Prediction");
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_ZMP_Bipedal_StateJump");
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_ZMP_Bipedal");

    // Com states
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_Com");
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_Com_Acc");

    // Com-vel states
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_Com_Vel_Prediction");
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_Com_Vel_StateJump");
    getControllerWithLogger_()->logger().removeLogEntry("FloatingBaseState_Com_Vel");
  }
  else
  {
     std::runtime_error("Trying to remove loggers without adding them!");
  }
}

void ImpactAwareFloatingBaseConstraint::addMcContactGuiItems(mc_control::fsm::Controller & ctl) const
{
  contactSetPtr_->addGuiItems(ctl);
}

void ImpactAwareFloatingBaseConstraint::addMcAreasGuiItems(mc_control::fsm::Controller & ctl) const
{
  switch(getConstrainingStatus())
  {
    case 1:
      if(getParams().updateMcZMPArea)
      {
        getMcZMPArea()->addGuiItems(ctl);
      }
      break;
    case 3:
      if(getParams().updateMcDCMArea)
      {
        getMcZMPArea()->addGuiItems(ctl);
        getMcComArea()->addGuiItems(ctl);
        getMcDCMArea()->addGuiItems(ctl);
      }
      break;
    case 4:

      if(getParams().updateMcZMPArea)
      {
        getMcZMPArea()->addGuiItems(ctl);
      }
      if(getParams().updateMcDCMArea)
      {
        getMcComArea()->addGuiItems(ctl);
        getMcDCMArea()->addGuiItems(ctl);
      }
      break;
    default:
      throw std::runtime_error("The constraining status is not set correctly");
  }
}

void ImpactAwareFloatingBaseConstraint::addFloatingBaseGuiItems(mc_control::fsm::Controller & ctl) const
{
  mc_rtc::gui::ArrowConfig arrow_config_zmp({0., 1., 0.});
  arrow_config_zmp.start_point_scale = 0.0;
  arrow_config_zmp.end_point_scale = 0.0;

  mc_rtc::gui::ArrowConfig arrow_config_dcm({1., 0., 0.});
  arrow_config_dcm.start_point_scale = 0.0;
  arrow_config_dcm.end_point_scale = 0.0;
  arrow_config_dcm.shaft_diam *= 2.0;
  arrow_config_dcm.head_diam *= 1.5;

  ctl.gui()->addElement({"DCM"}, mc_rtc::gui::Arrow("DCM", arrow_config_dcm,
                                                    [this]() {
                                                      const auto & c = robot().com();
                                                      const auto & cdot = robot().comVelocity();
                                                      double w = std::sqrt(9.81 / 0.78);
                                                      Eigen::Vector3d dcm3 = (c + cdot / w); // dcmImpulse_->dcm;
                                                      return Eigen::Vector3d{dcm3.x(), dcm3.y(), 0.};
                                                    },
                                                    [this]() {
                                                      // End of the arrow
                                                      return robot().com();
                                                    }));

  ctl.gui()->addElement({"ZMP"}, mc_rtc::gui::Arrow("Arrow", arrow_config_zmp,
                                                    [this]() { return getFloatingBaseStates().mcZMP.current; },
                                                    [this]() {
                                                      // End of the arrow
                                                      return robot().com();
                                                    }));
}

} // end of namespace mc_impact
