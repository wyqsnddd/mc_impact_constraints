#include "ImpactAwareCOMVelConstraint.h"


namespace mc_impact{
ImpactAwareCOMVelConstraint::ImpactAwareCOMVelConstraint(std::shared_ptr<mi_qpEstimator> predictorPtr,
		const mc_rbdyn::Robot & realRobot,
		const ImpactAwareConstraintParams<Eigen::Vector2d> & params)
	:mc_solver::InequalityConstraintRobot(predictorPtr->getSimRobot().robotIndex()), predictorPtr_(predictorPtr), params_(params), robot_(predictorPtr->getSimRobot()), realRobot_(realRobot)
{
  // Update the constraint size
  A_.resize(nrInEq(), dof_());
  A_.setZero();
  b_.resize(nrInEq());
  b_.setZero();

  // Initialize the com Jacobian
  comJacobianPtr_ = std::make_shared<rbd::CoMJacobian>(robot_.mb());

  // Initialize the support polygon
  
  poleOne_ = -1.0;
  poleTwo_ = -0.2;
  updateGains_();

  std::cout << red << "Created ImpactAwareCOMVelConstraint." << reset << std::endl;
}// end of constructor


void ImpactAwareCOMVelConstraint::compute()
{

  //robotJointVelocity_ = (rbd::dofToVector(realRobot().mb(), realRobot().mbc().alpha));
  robotJointVelocity_ = (rbd::dofToVector(robot().mb(), robot().mbc().alpha));

  // update the sup 
  
    // Update the COM vel constraint:
  updateCOMVelConstraint_();
}

void ImpactAwareCOMVelConstraint::updateCOMVelBounds_()
{

  updateGains_();
  // Here we assume that com_ref and com_dot_ref are all zeros.
  COMStates_.ComVelXUpperBound = (zUpperBound_ - COMStates_.Com.x()*pGain_)/dGain_; 
  COMStates_.ComVelXLowerBound = (zLowerBound_ - COMStates_.Com.x()*pGain_)/dGain_; 

  if(getParams().debug)
  {

    std::cout<<"The COM is: "<<red<<  COMStates_.Com.transpose() << reset<<", and the omege is: "<<red<<getOmega() <<reset<<std::endl;
    std::cout<<"The pGain is: "<<red<<  pGain_ << reset<<", and the dGain is: "<<red<<dGain_<<reset<<std::endl;
    std::cout<<"The zmp upper bound is: "<<red<<  zUpperBound_ << reset<<", and the zmp lower bound is: "<<red<<zLowerBound_<<reset<<std::endl;
    std::cout<<"The upper bound is: "<<red<<  COMStates_.ComVelXUpperBound << reset<<", and the lower is: "<<red<<COMStates_.ComVelXLowerBound <<reset<<std::endl;
  }
}

void ImpactAwareCOMVelConstraint::updateCOMVelConstraint_()
{
  // Update omega for DCM calculation
  calcOmega_(robot().com().z());

  Eigen::MatrixXd comJacobian = comJacobianPtr_->jacobian(robot().mb(), robot().mbc());
  Eigen::MatrixXd comJacobianDot = comJacobianPtr_->jacobianDot(robot().mb(), robot().mbc());

  Eigen::VectorXd alphaD = rbd::dofToVector(robot().mb(), robot().mbc().alphaD);

  COMStates_.Com = realRobot().com();

  //COMStates_.ComAcc = realRobot().comAcceleration();
  // The realRobot acceleration is not usable.
  
  // COMStates_.ComAcc = robot().comAcceleration();
  // The simulated robot COM acceleration is extremely small.
  
  COMStates_.ComAcc = comJacobian*alphaD + comJacobianDot*robotJointVelocity_;

  // COM velocity

  //COMStates_.ComVel.current = realRobot().comVelocity();
  COMStates_.ComVel.current = robot().comVelocity();
  if(getParams().debug)
  {
    std::cout << red << "ImpactAwareCOMVelConstraint computed FloatingBase-States." << reset << std::endl;
  }
  const Eigen::MatrixXd &  cmmMatrix = getPredictor()->getCmm()->matrix();

  double mass_inv = 1.0/robot().mass();

  // Special Jac one
  
  // Note that the CMM matrix starts with the angular part and then translation, we need to shift the row number for 3.
  Eigen::VectorXd specialJacOne_cmm = mass_inv * cmmMatrix.block(0, 0, 1, dof_()) * getPredictor()->getJacobianDeltaAlpha();
  Eigen::VectorXd specialJacOne =  comJacobian.block(0, 0, 1, dof_()) * getPredictor()->getJacobianDeltaAlpha();

  Eigen::VectorXd specialJacTwo_cmm = mass_inv * cmmMatrix.block(0, 0, 1, dof_()) * getPredictor()->getJacobianTwoDeltaAlpha();
  Eigen::VectorXd specialJacTwo = comJacobian.block(0, 0, 1, dof_()) * getPredictor()->getJacobianTwoDeltaAlpha();


  
  // Special Jac Two
  // Saggital(X)-direction: 
  A_.block(0, 0, 1, dof_()) = getParams().dt * specialJacOne_cmm; 
  A_.block(1, 0, 1, dof_()) = - A_.block(0, 0, 1, dof_());
  //A_.block(1, 0, 1, dof_()) = - getParams().dt * specialJacOne; 

  // Update the COM vel bounds:
  updateCOMVelBounds_();

  if(getParams().debug)
  {
    Eigen::Matrix3Xd specialJacOne_cmm_full = mass_inv * cmmMatrix.block(0, 0, 3, dof_()) * getPredictor()->getJacobianDeltaAlpha();

    Eigen::Matrix3Xd specialJacTwo_cmm_full = mass_inv * cmmMatrix.block(0, 0, 3, dof_()) * getPredictor()->getJacobianTwoDeltaAlpha();
  

  Eigen::VectorXd specialJacOne_full =  comJacobian * getPredictor()->getJacobianDeltaAlpha();

  Eigen::VectorXd specialJacTwo_full  = comJacobian * getPredictor()->getJacobianTwoDeltaAlpha();

    temp_comVel_jump_comJacobian_acc = getParams().dt*specialJacOne_full*alphaD + specialJacTwo_full*robotJointVelocity_;

    temp_comVel_jump_comJacobian = comJacobian*getPredictor()->getJointVelJump();

    temp_comVel_jump_cmm = mass_inv * cmmMatrix.block(0, 0, 3, dof_())*getPredictor()->getJointVelJump();

    temp_comVel_jump_cmm_two = mass_inv * cmmMatrix.block(3, 0, 3, dof_())*getPredictor()->getJointVelJump();

    COMStates_.ComVel.stateJump =  
	 getParams().dt*specialJacOne_cmm_full*alphaD + specialJacTwo_cmm_full.transpose() * robotJointVelocity_;

    COMStates_.ComVel.oneStepPreview = COMStates_.ComVel.current + COMStates_.ComVel.stateJump;


    std::cout << red << "COM Vel upper bound is: " << getCOMStates().ComVelXUpperBound<< reset << std::endl;
    std::cout << red << "currnt COM vel is: " << COMStates_.ComVel.current.x()<< reset << std::endl;
    std::cout << red << "currnt COM acc is: " << COMStates_.ComAcc.x()<< reset << std::endl;
    std::cout << red << "currnt dt is: " << getParams().dt << reset << std::endl;

    Eigen::VectorXd alphaD = rbd::dofToVector(robot().mb(), robot().mbc().alphaD);

    std::cout << red << "lhs is: " << A_.block(0, 0, 1, dof_())*alphaD << reset << std::endl;
    std::cout << red << "rhs is: " << - specialJacTwo.transpose() * robotJointVelocity_ << reset << std::endl;
    

    double com_vel_jump = (getParams().dt*specialJacOne*alphaD)(0) + specialJacTwo.transpose() * robotJointVelocity_;
    //double com_vel_jump_second = (getParams().dt*specialJacOne*alphaD)(0) + specialJacOne.transpose() * robotJointVelocity_;

    double com_vel_jump_cmm = (getParams().dt*specialJacOne_cmm*alphaD)(0) + specialJacTwo_cmm.transpose() * robotJointVelocity_;
    //double com_vel_jump_cmm_two = (getParams().dt*specialJacOne_cmm*alphaD)(0) + specialJacOne_cmm.transpose() * robotJointVelocity_;

    double com_vel_next = COMStates_.ComVel.current.x() + COMStates_.ComAcc.x()*getParams().dt;

    std::cout << yellow << "The predicted next COM-vel is: " << com_vel_next << reset << std::endl;

    std::cout << yellow << "The predicted COM-vel jump is: " << com_vel_jump << reset << std::endl;
    std::cout << yellow << "The predicted post-impact COM vel is: " << com_vel_jump + com_vel_next  << reset << std::endl;

    //std::cout << yellow << "The second predicted COM-vel jump is: " << com_vel_jump_second << reset << std::endl;
    //std::cout << yellow << "The second predicted post-impact COM vel is: " << com_vel_jump_second + com_vel_next  << reset << std::endl;

    std::cout << cyan << "The cmm predicted COM-vel jump is: " << com_vel_jump_cmm << reset << std::endl;
    std::cout << cyan << "The cmm predicted post-impact COM vel is: " << com_vel_jump_cmm + com_vel_next  << reset << std::endl;

    //std::cout << cyan << "The cmm second predicted COM-vel jump is: " << com_vel_jump_cmm_two<< reset << std::endl;
    //std::cout << cyan << "The cmm second predicted post-impact COM vel is: " << com_vel_jump_cmm_two+ com_vel_next  << reset << std::endl;

  }

  b_(0) = getCOMStates().ComVelXUpperBound  
	  - COMStates_.ComVel.current.x()
	  - COMStates_.ComAcc.x()*getParams().dt
	  - specialJacTwo_cmm.transpose() * robotJointVelocity_; 
  //b_(1) = - getCOMStates().ComVelXLowerBound - (b_(0) - getCOMStates().ComVelXUpperBound);
  b_(1) =  - getCOMStates().ComVelXLowerBound  
	  + COMStates_.ComVel.current.x()
	  + COMStates_.ComAcc.x()*getParams().dt
	  + specialJacTwo_cmm.transpose() * robotJointVelocity_; 



  // Lateral(Y)-direction:
  
}


void ImpactAwareCOMVelConstraint::logConstraint(mc_control::fsm::Controller & ctl) 
{

  setControllerWithLogger_(ctl);
  ctl.logger().addLogEntry("ComVel_realrobot", [this]() {
    return realRobot().comVelocity(); 
  });

  setControllerWithLogger_(ctl);
  ctl.logger().addLogEntry("ComVel_x_upper_bound", [this]() {
    return getCOMStates().ComVelXUpperBound;
  });

  ctl.logger().addLogEntry("ComVel_x_lower_bound", [this]() {
    return getCOMStates().ComVelXLowerBound;
  });
  
  ctl.logger().addLogEntry("ComVel_y_upper_bound", [this]() {
    return getCOMStates().ComVelYUpperBound;
  });

  ctl.logger().addLogEntry("ComVel_y_lower_bound", [this]() {
    return getCOMStates().ComVelYLowerBound;
  });


 if(getParams().debug)
  { 

  ctl.logger().addLogEntry("ComVel_predicted_jump_comJacobian_acc", [this]() {
    return temp_comVel_jump_comJacobian_acc; 
  });


  ctl.logger().addLogEntry("ComVel_predicted_jump_comJacobian", [this]() {
    return temp_comVel_jump_comJacobian; 
  });

  ctl.logger().addLogEntry("ComVel_predicted_jump_cmm", [this]() {
    return temp_comVel_jump_cmm ;
  });

  ctl.logger().addLogEntry("ComVel_predicted_jump", [this]() {
    return getCOMStates().ComVel.stateJump;
  });

  ctl.logger().addLogEntry("ComVel_next_preImpact", [this]() {
    return  (Eigen::Vector3d)(COMStates_.ComVel.current + COMStates_.ComAcc*getParams().dt);
  });


  }

  setLogger_();
}

void ImpactAwareCOMVelConstraint::removeLog() 
{

  if(checkLoggersAdded_())
  {
    getControllerWithLogger_()->logger().removeLogEntry("ComVel_x_upper_bound");
    getControllerWithLogger_()->logger().removeLogEntry("ComVel_x_lower_bound");
    getControllerWithLogger_()->logger().removeLogEntry("ComVel_y_upper_bound");
    getControllerWithLogger_()->logger().removeLogEntry("ComVel_y_lower_bound");

    if(getParams().debug)
    {
      getControllerWithLogger_()->logger().removeLogEntry("COMVel_predicted_jump");
      getControllerWithLogger_()->logger().removeLogEntry("ComVel_predicted_jump_comJacobian");
      getControllerWithLogger_()->logger().removeLogEntry("ComVel_predicted_jump_cmm");
      getControllerWithLogger_()->logger().removeLogEntry("ComVel_predicted_jump_comJacobian_acc");

    }
  }
  else
  {
     std::runtime_error("Trying to remove loggers without adding them!");
  }
}



} // namespace mc_impact
