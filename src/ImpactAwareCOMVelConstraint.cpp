#include "ImpactAwareCOMVelConstraint.h"


namespace mc_impact{
ImpactAwareCOMVelConstraint::ImpactAwareCOMVelConstraint(std::shared_ptr<mi_qpEstimator> predictorPtr,
		const ImpactAwareConstraintParams<Eigen::Vector2d> & params)
	:mc_solver::InequalityConstraintRobot(predictorPtr->getSimRobot().robotIndex()), predictorPtr_(predictorPtr), params_(params), robot_(predictorPtr->getSimRobot())
{
  // Update the constraint size
  A_.resize(nrInEq(), dof_());
  A_.setZero();
  b_.resize(nrInEq());
  b_.setZero();

  // Initialize the com Jacobian
  comJacobianPtr_ = std::make_shared<rbd::CoMJacobian>(robot_.mb());

  // Initialize the support polygon
  
  poleOne_ = -1;
  poleTwo_ = -0.2;
  updateGains_();

  std::cout << red << "Created ImpactAwareCOMVelConstraint." << reset << std::endl;
}// end of constructor


void ImpactAwareCOMVelConstraint::compute()
{

  robotJointVelocity_ = (rbd::dofToVector(robot().mb(), robot().mbc().alpha));

  // update the sup 
  
  // Update the COM vel bounds:
  updateCOMVelBounds_();

  // Update the COM vel constraint:
  updateCOMVelConstraint_();
}

void ImpactAwareCOMVelConstraint::updateCOMVelBounds_()
{

  updateGains_();
  // Here we assume that com_ref and com_dot_ref are all zeros.
  COMStates_.ComVelXUpperBound = (zUpperBound_ - COMStates_.Com.x()*pGain_)/dGain_; 
  COMStates_.ComVelXLowerBound = (zLowerBound_ - COMStates_.Com.x()*pGain_)/dGain_; 
}

void ImpactAwareCOMVelConstraint::updateCOMVelConstraint_()
{
  // Update omega for DCM calculation
  calcOmega_(robot().com().z());

  Eigen::MatrixXd comJacobian = comJacobianPtr_->jacobian(robot().mb(), robot().mbc());

  COMStates_.Com = robot().com();
  COMStates_.ComAcc = robot().comAcceleration();

  // COM velocity

  COMStates_.ComVel.current = robot().comVelocity();
  COMStates_.ComVel.stateJump = comJacobian * getPredictor()->getJointVelJump();
  COMStates_.ComVel.oneStepPreview = COMStates_.ComVel.current + COMStates_.ComVel.stateJump;

  if(getParams().debug)
  {
    std::cout << red << "ImpactAwareCOMVelConstraint computed FloatingBase-States." << reset << std::endl;
  }
  const Eigen::MatrixXd &  cmmMatrix = getPredictor()->getCmm()->matrix();

  double mass_inv = 1.0/robot().mass();

  // Special Jac one
  
  Eigen::VectorXd specialJacOne = mass_inv * cmmMatrix.block(0, 0, 1, dof_()) * getPredictor()->getJacobianDeltaAlpha();

  Eigen::VectorXd specialJacTwo = mass_inv * cmmMatrix.block(0, 0, 1, dof_()) * getPredictor()->getJacobianTwoDeltaAlpha();

  // Special Jac Two
  // Saggital(X)-direction: 
  A_.block(0, 0, 1, dof_()) = getParams().dt * specialJacOne; 
  A_.block(1, 0, 1, dof_()) = - A_.block(0, 0, 1, dof_());

  b_(0) = getCOMStates().ComVelXUpperBound
	  - COMStates_.ComVel.oneStepPreview.x()
	  - specialJacTwo.transpose() * robotJointVelocity_; 
  b_(1) = getCOMStates().ComVelXLowerBound - (b_(0) - getCOMStates().ComVelXUpperBound);

  // Lateral(Y)-direction:
  
}


} // namespace mc_impact
