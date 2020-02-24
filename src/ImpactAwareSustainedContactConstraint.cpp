#include "ImpactAwareSustainedContactConstraint.h"

namespace mc_impact
{

ImpactAwareSustainedContact::ImpactAwareSustainedContact(
                               mi_qpEstimator & predictor,
		 const McContactParams & mcContactParams,
		 const ImpactAwareConstraintParams<Eigen::Vector2d> & constraintParams)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), 
  predictor_(predictor), constraintParams_(constraintParams), mcContactParams_(mcContactParams)
{

  initializeImpactAwareCopConstraint_();
  initializeImpactAwareFrictionConstraint_();

  measuredWrench_.couple().setZero();
  measuredWrench_.force().setZero();

  robotJointVelocity_.resize(dof());

  A_.resize(nrInEq(), dof());
  A_.setZero();
  b_.resize(nrInEq());
  b_.setZero();
}

void ImpactAwareSustainedContact::initializeImpactAwareFrictionConstraint_()
{
  multiplier_.resize(2, 3);
  multiplier_.setZero();
  multiplier_(0, 0) = 1.0;
  multiplier_(0, 2) = -getParams().frictionCoe;
  multiplier_(1, 1) = 1.0;
  multiplier_(1, 2) = -getParams().frictionCoe;

  A_friction_.resize(frictionConstraintSize_(), dof());
  b_friction_.resize(frictionConstraintSize_());

  A_friction_.setZero();
  b_friction_.setZero();
}
void ImpactAwareSustainedContact::initializeImpactAwareCopConstraint_()
{
  G_cop_ = Eigen::MatrixXd::Zero(4, 6);

  double maxX = getParams().halfX;
  double minX = -getParams().halfX;
  double maxY = getParams().halfY;
  double minY = -getParams().halfY;

  // - cy - max_x * fz
  G_cop_(0, 1) = -1;
  G_cop_(0, 5) = -maxX;

  // cy + min_x * fz
  G_cop_(1, 1) = 1;
  G_cop_(1, 5) = minX;

  // cx - max_y * fz
  G_cop_(2, 0) = 1;
  G_cop_(2, 5) = -maxY;

  // - cx + min_y * fz
  G_cop_(3, 0) = -1;
  G_cop_(3, 5) = minY;

  A_cop_.resize(copConstraintSize_(), dof());
  A_cop_.setZero();

  b_cop_.resize(copConstraintSize_());
  b_cop_.setZero();

}

void ImpactAwareSustainedContact::updateCopLogs_()
{

  // (1) Update the COP 
  // In the body frame, the contact surface normal is always the Z axis.
  cop_ = copCalculation(Eigen::Vector3d::UnitZ(), measuredWrench());

  // (2) Update the perturbed COP
   
  
  auto perturbedWrench = measuredWrench();
  perturbedWrench.force() += predictor_.getEndeffector(getParams().bodyName).estimatedAverageImpulsiveForce;

  cop_stateJump_ = copCalculation(Eigen::Vector3d::UnitZ(), perturbedWrench);

  /*
  cop_perturb_whole_.setZero();

  double perturbedNormalForce_whole =
      perturbedNormalForce + predictor_.getEndeffector(getParams().bodyName).perturbedWrench.force().z();

  cop_perturb_whole_.x() =
      -(measuredWrench().couple().y() + predictor_.getEndeffector(getParams().bodyName).perturbedWrench.couple().y())
      / perturbedNormalForce_whole;

  cop_perturb_whole_.y() =
      (measuredWrench().couple().x() + predictor_.getEndeffector(getParams().bodyName).perturbedWrench.couple().x())
      / perturbedNormalForce_whole;
      */

}
void ImpactAwareSustainedContact::computeImpactAwareCopConstraint_()
{
  const auto & J_deltaF = predictor_.getJacobianDeltaF(getParams().bodyName);

  A_cop_ = (getConstraintParams().dt / getConstraintParams().impactDuration) * G_cop_.block(0, 3, 4, 3) * J_deltaF;


  b_cop_ = -(G_cop_ * measuredWrench().vector()
         + G_cop_.block(0, 3, 4, 3) * J_deltaF * robotJointVelocity_/ getConstraintParams().impactDuration);

}

void ImpactAwareSustainedContact::computeImpactAwareFrictionConstraint_()
{
  const auto & J_deltaF = predictor_.getJacobianDeltaF(getParams().bodyName);

  A_friction_ = (getConstraintParams().dt / getConstraintParams().impactDuration) * multiplier_ * J_deltaF;

  b_friction_ = -multiplier_
       * (predictor_.getSimRobot().bodyWrench(getParams().bodyName).force()
          + J_deltaF * robotJointVelocity_/ getConstraintParams().impactDuration);

}
void ImpactAwareSustainedContact::compute()
{

   // Update the measured Wrench 
    measuredWrench_ = robot().bodyWrench(getParams().bodyName);

   // Update robot joint velocity: 
   rbd::paramToVector(robot().mbc().alpha, robotJointVelocity_);

   // Update cop 
   computeImpactAwareCopConstraint_();
   updateCopLogs_();

   // Update friction
   computeImpactAwareFrictionConstraint_();

   // Evaluate the matrix A_ and vector b_;
   A_.block(0, 0, copConstraintSize_(), dof()) = A_cop_;
   A_.block(copConstraintSize_(), 0, frictionConstraintSize_(), dof()) = A_friction_;

   b_.segment(0, copConstraintSize_()) = b_cop_; 
   b_.segment(copConstraintSize_(), frictionConstraintSize_()) = b_friction_; 

}


Eigen::Matrix3d ImpactAwareSustainedContact::crossProduct(const Eigen::Vector3d & input)
{

  Eigen::Matrix3d skewSymmetricMatrix = Eigen::Matrix3d::Zero();

  skewSymmetricMatrix(0, 1) = -input(2);
  skewSymmetricMatrix(1, 0) = input(2);

  skewSymmetricMatrix(0, 2) = input(1);
  skewSymmetricMatrix(2, 0) = -input(1);

  skewSymmetricMatrix(1, 2) = -input(0);
  skewSymmetricMatrix(2, 1) = input(0);

  return skewSymmetricMatrix;
}

Eigen::Vector3d ImpactAwareSustainedContact::copCalculation(const Eigen::Vector3d & surfaceNormal, const sva::ForceVecd inputWrench)
{

  Eigen::Matrix3d crossGroundNormal = mc_impact::crossMatrix(surfaceNormal);

  double denominator = surfaceNormal.dot(inputWrench.force());

  return (crossGroundNormal * inputWrench.couple()) * (1.0 / denominator);
}


} // namespace mc_impact
