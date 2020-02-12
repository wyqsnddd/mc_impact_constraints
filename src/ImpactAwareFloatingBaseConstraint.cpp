# include "ImpactAwareFloatingBaseContraint.h"


namespace mc_impact
{
   ImpactAwareFloatingBaseConstraint::ImpactAwareFloatingBaseConstraint(mi_qpEstimator & predictor,
                 const mc_rbdyn::Robot & realRobot,
                 std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr,
                 std::shared_ptr<mc_impact::McComArea> mcComAreaPtr,
                 const ImpactAwareConstraintParams<Eigen::Vector2d> & params)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), predictor_(predictor), realRobot_(realRobot),
  mcZMPAreaPtr_(mcZMPAreaPtr), mcComAreaPtr_(mcComAreaPtr), params_(params)
{

  mcDCMAreaPtr_ = std::make_shared<mc_impact::McDCMArea>(mcZMPAreaPtr_, mcComAreaPtr_);


  // Set the constraining status of the floating-base state. 
  constrainingStatus_ = 0;
  if (getParams().constrainingZMP)
  {
    constrainingStatus_ +=1; 
  }
  if (getParams().constrainingZMP)
  {
    constrainingStatus_ +=3; 
  }

  // Initialize the com Jacobian 
  comJacobianPtr_ = std::make_shared<rbd::CoMJacobian>(realRobot_.mb());

}
void ImpactAwareFloatingBaseConstraint::reset_()
{
  int numVertexZMP = mcZMPAreaPtr_->getNumVertex();

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


  //int numVertexDCM = mcDCMAreaPtr_->getNumVertex();
  //A_dcm_ = Eigen::MatrixXd::Zero(numVertexDCM, 6);
  // Modify the ieqConstraintBlocks.
  pointsToInequalityMatrix(getParams().dcmAreaVertexSet, ieqConstraintBlocksDCM_.G, ieqConstraintBlocksDCM_.h,
                           getParams().lowerSlope, getParams().upperSlope);
  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  //A_dcm_.block(0, 0, numVertexDCM, 1) = getIeqBlocksDCM().G.block(0, 1, numVertexDCM, 1);
  /// A(:,1) = -G_x
  //A_dcm_.block(0, 1, numVertexDCM, 1) = -getIeqBlocksDCM().G.block(0, 0, numVertexDCM, 1);
  /// A(:,5) = h
  //A_dcm_.block(0, 5, numVertexDCM, 1) = -getIeqBlocksDCM().h;


  // Initialie the building blocks
  robotJointVelocity_.resize(dof_());
  A_.resize(nrInEq(), dof_());
  b_.resize(nrInEq());

}
void ImpactAwareFloatingBaseConstraint::updateMcZMPAreas_(double height)
{

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

void ImpactAwareFloatingBaseConstraint::updateMcDCMAreas_()
{
  // --------------- (2) Update the Multi-contact ZMP area:
  
  // If ZMP constraint is not used, we also need to update the ZMP area: 
  if(zmpConstraintEnabled_())
  {
    mcZMPAreaPtr_->updateMcZMPArea(2.0);
  }

  mcComAreaPtr_->updateMcComArea();
  mcDCMAreaPtr_->updateMcDCMArea();
  // int numVertex = static_cast<int>(iniVertexSet_.size());
  //int numVertexDCM = getMcDCMArea()->getNumVertex();

  //A_dcm_ = Eigen::MatrixXd::Zero(numVertexDCM, 6);

  // Set the inequality matrix blocks
  setIeqBlocksDCM(getMcDCMArea()->getIeqConstraint());

  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  //A_dcm_.block(0, 0, numVertexDCM, 1) = getIeqBlocksDCM().G.block(0, 1, numVertexDCM, 1);
  /// A(:,1) = -G_x
  //A_dcm_.block(0, 1, numVertexDCM, 1) = -getIeqBlocksDCM().G.block(0, 0, numVertexDCM, 1);
  /// A(:,5) = h
  //A_dcm_.block(0, 5, numVertexDCM, 1) = -getIeqBlocksDCM().h;

}

void ImpactAwareFloatingBaseConstraint::getZMPBlocks(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench)
{
  exWrench.setZero();
  //int dof = predictor_.getSimRobot().mb().nrDof();
  //int dof = realRobot().mb().nrDof();
  sumJac.resize(6, dof_());
  sumJac.setZero();
  // sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  // (1) Go through the bodies with contact

  // for(auto idx = getParams().contactSetPtr->getContactMap().begin(); idx != getParams().contacts.end(); ++idx)

  for(auto & contactPair : getParams().contactSetPtr->getContactMap())
  {
    std::string bodyName = contactPair.second.getContactParams().bodyName;
    //sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(contactPair.second.getContactParams().bodyName).inv();
    sva::PTransformd X_ee_0 = realRobot().bodyPosW(contactPair.second.getContactParams().bodyName).inv();
    // sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(idx->bodyName).inv();
    // exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().forceSensor(idx->sensorName).wrench().vector();
    //exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(bodyName).vector();
    exWrench += X_ee_0.dualMatrix() * realRobot().bodyWrench(bodyName).vector();

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
      //sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();
      sva::PTransformd X_ee_0 = realRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();

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

void ImpactAwareFloatingBaseConstraint::calculateZMP_()
{

  // We use the joint velocity used by the impact model
  Eigen::VectorXd jointVel = predictor_.getImpactModels().begin()->second->getJointVel();
  double inv_t = (1 / predictor_.getImpactModels().begin()->second->getImpactDuration());

  //Eigen::Vector6d local_exWrench;
  //local_exWrench.setZero();

  Eigen::Vector6d contactWrench;
  //Eigen::Vector6d impactWrench;
  contactWrench.setZero();
  //impactWrench.setZero();

  Eigen::MatrixXd contactWrenchJac; // Forces due to established contacts.
  Eigen::MatrixXd  impactWrenchJac; // Predicted forces due to expected impacts.
   
  contactWrenchJac.resize(6, dof_());
  contactWrenchJac.setZero();

  impactWrenchJac.resize(6, dof_());
  impactWrenchJac.setZero();


  // (1) Go through all the contact
  // Collect all the local Jacobians from all the contacts 
  for(auto & contactPair : getParams().contactSetPtr->getContactMap())
  {

    std::string bodyName = contactPair.second.getContactParams().bodyName;
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(bodyName).inv();
    if(contactPair.second.inContact()){
      // contact established
      contactWrench += X_ee_0.dualMatrix() * realRobot().bodyWrench(bodyName).vector();
      contactWrenchJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(bodyName);
    }else{
       // contact not established, impact is expected.
      //impactWrench += X_ee_0.dualMatrix() * realRobot().bodyWrench(bodyName).vector();
      impactWrenchJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(bodyName);
    }
  } // end of for

  // Sum of the impulse forces at the established contacts:
  Eigen::Vector6d impulseContactWrench;
  impulseContactWrench.setZero();
  impulseContactWrench= inv_t * contactWrenchJac* jointVel;
  //double denominator = contactWrench(5) + impulseWrenchSum(5);
  Eigen::Vector6d contactWrenchSum = impulseContactWrench + contactWrench;

  sva::ForceVecd sva_contactWrenchSum;
  sva_contactWrenchSum.force() = contactWrenchSum.segment(3, 3); // Force starts at 3.
  sva_contactWrenchSum.couple() = contactWrenchSum.segment(0, 3); // Torque starts at 0.

  floatingBaseStates_.bipedalZMP.oneStepPreview = getMcZMPArea()->zmpCalculation(Eigen::Vector3d::UnitZ(), sva_contactWrenchSum); 
  floatingBaseStates_.bipedalZMP.current = getMcZMPArea()->getBipedalZMP();
  floatingBaseStates_.bipedalZMP.stateJump =  floatingBaseStates_.bipedalZMP.oneStepPreview - floatingBaseStates_.bipedalZMP.current; 


  Eigen::Vector6d sumWrench;
  sumWrench.setZero();
  // Add the predicted impulse 
  sumWrench = inv_t * impactWrenchJac* jointVel + contactWrenchSum;

  sva::ForceVecd sva_sum;
  sva_sum.force() = sumWrench.segment(3, 3); // Force starts at 3.
  sva_sum.couple() = sumWrench.segment(0, 3); // Torque starts at 0.

  floatingBaseStates_.mcZMP.oneStepPreview = getMcZMPArea()->zmpCalculation(Eigen::Vector3d::UnitZ(), sva_sum); 
  floatingBaseStates_.mcZMP.current = getMcZMPArea()->getMcZMP();
  floatingBaseStates_.mcZMP.current = floatingBaseStates_.mcZMP.oneStepPreview - floatingBaseStates_.mcZMP.current; 

}

void ImpactAwareFloatingBaseConstraint::updateFloatingBaseState_()
{
  // Update omega for DCM calculation
  calcOmega(realRobot().com().z());

  Eigen::MatrixXd comJacobian = comJacobianPtr_->jacobian(realRobot().mb(), realRobot().mbc());
  Eigen::MatrixXd dcmJacobian  = (comJacobian * predictor_.getJacobianDeltaAlpha()).block(0, 0, 2, dof_());

  floatingBaseStates_.Com = realRobot().com();
  floatingBaseStates_.ComAcc = realRobot().comAcceleration();

  // (1) COM velocity
  
  floatingBaseStates_.ComVel.current = realRobot().comVelocity();
  floatingBaseStates_.ComVel.stateJump = comJacobian * predictor_.getJointVelJump();
  floatingBaseStates_.ComVel.oneStepPreview = floatingBaseStates_.ComVel.current + floatingBaseStates_.ComVel.stateJump;

  // (2) ZMP 
  calculateZMP_();

  // (3) DCM
  floatingBaseStates_.DCM.current.segment(0,2) = floatingBaseStates_.Com.segment(0, 2) + floatingBaseStates_.ComVel.current.segment(0, 2) / getOmega();
  floatingBaseStates_.DCM.stateJump.segment(0,2) = comJacobian.block(0, 0, 2, dof_()) * predictor_.getJointVelJump() / getOmega();
  floatingBaseStates_.DCM.oneStepPreview.segment(0,2) =  floatingBaseStates_.DCM.current.segment(0,2) +  floatingBaseStates_.DCM.stateJump.segment(0,2);


}

void ImpactAwareFloatingBaseConstraint::updateDCMConstraint_()
{
  // Calculates the new DCM constraint blocks internally.  
  updateMcDCMAreas_();

   // Should we use the real robot or the simulated one? 
  //const auto & robot = predictor_.getSimRobot();
  //int dof = dof_(); 

  //Eigen::MatrixXd jacCom= comJacobianPtr_->jacobian(realRobot().mb(), realRobot().mbc());
  // std::cout<<"comJacobian size is: "<<comJacobian.rows() << ", "<<comJacobian.cols()<<std::endl;
  
  Eigen::MatrixXd jacDCM = (comJacobianPtr_->jacobian(realRobot().mb(), realRobot().mbc())* predictor_.getJacobianDeltaAlpha()).block(0, 0, 2, dof_());

  A_.block(getDCMRowNr_(), 0, getMcDCMArea()->getNumVertex(), 6) = getParams().dt / getOmega() * ieqConstraintBlocksDCM_.G * jacDCM;

  b_.segment(getDCMRowNr_(), getMcDCMArea()->getNumVertex()) = (ieqConstraintBlocksDCM_.h - ieqConstraintBlocksDCM_.G * floatingBaseStates_.DCM.current) - ieqConstraintBlocksDCM_.G * jacDCM * robotJointVelocity_/ getOmega();

}
void ImpactAwareFloatingBaseConstraint::updateZMPConstraint_()
{
  Eigen::MatrixXd sumJacZMP;
  Eigen::Vector6d sumWrenchZMP;
  getZMPBlocks(sumJacZMP, sumWrenchZMP);
  updateMcZMPAreas_(2.0);

  // We assume that we always start from ZMP constraint unless it is not used. 
  A_.block(0, 0, getMcZMPArea()->getNumVertex(), 6) = (getParams().dt / getParams().impactDuration) * A_zmp_ * sumJacZMP;

  // Read the robot joint velocities.
  //rbd::paramToVector(robot.mbc().alpha, robotJointVelocity_);
  b_.segment(0, getMcZMPArea()->getNumVertex()) = -(A_zmp_ * sumWrenchZMP + A_zmp_ * sumJacZMP * robotJointVelocity_/ getParams().dt);
}

void ImpactAwareFloatingBaseConstraint::compute()
{
   // Update the ZMP, ComVel and DCM states
   updateFloatingBaseState_();

   //const auto & robot = predictor_.getSimRobot();
   // Read the robot joint velocity:
   rbd::paramToVector(realRobot().mbc().alpha, robotJointVelocity_);
   switch(getConstrainingStatus()){
   case 1:
      updateZMPConstraint_();
   case 3:
      updateDCMConstraint_();
   case 4:
      updateZMPConstraint_();
      updateDCMConstraint_();
   default:
      throw std::runtime_error("The constraining status is not set correctly");
   }
}



}



