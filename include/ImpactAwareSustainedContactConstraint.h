#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ConstraintUtils.h"

namespace mc_impact
{

struct ImpactAwareSustainedContact: public mc_solver::InequalityConstraintRobot
/*! \brief It contains two impact-aware constraints: (1) Center of pressure (2): friction within the friction cone.
 *
 */
{
  ImpactAwareSustainedContact(const std::string & mcContactName,
                 mi_qpEstimator & predictor,
                 const ImpactAwareConstraintParams<Eigen::Vector2d> & impactAwareConstraintParams);

  inline int maxInEq() const override
  {
    return nrInEq();
  }

  inline int nrInEq() const override
  {
    return copConstraintSize_() + frictionConstraintSize_();
  }

  inline std::string nameInEq() const override
  {
    return "ImpactAwareSustainedContact";
  }
  
  /*! \brief Measured COP in the body frame 
   */
  inline const Eigen::Vector3d & COP() const
  {
    return cop_;
  }

  /*! \brief Predicted state jump of the COP in the body frame 
   */
  inline const Eigen::Vector3d & CopStateJump() const
  {
    return cop_stateJump_;
  }

  /*
  inline const Eigen::Vector2d & getCopPerturbWhole() const
  {
    return cop_perturb_whole_;
  }
  */

  inline const Eigen::MatrixXd & A() const override
  {
    return A_;
  }
  inline const Eigen::VectorXd & bInEq() const override
  {
    return b_;
  }

  void compute() override;

  inline const std::string & getMcContactName() const
  {
    return mcContactName_;
  }

  inline const ImpactAwareConstraintParams<Eigen::Vector2d> & getConstraintParams() const
  {
    return constraintParams_;
  }

  inline const McContactParams & getParams() const
  {
    return constraintParams_.contactSetPtr->getContact(getMcContactName()).getContactParams();
  }

  const mc_rbdyn::Robot & robot() const
  {
    return predictor_.getSimRobot();
  }

  int dof() const
  {
    return  robot().mb().nrDof();
  }

  const sva::ForceVecd measuredWrench() const
  {
    return measuredWrench_; 
  }

  static Eigen::Vector3d copCalculation(const Eigen::Vector3d & normal, const sva::ForceVecd Wrench);

  static Eigen::Matrix3d crossProduct(const Eigen::Vector3d & input);


private:

  int copConstraintSize_() const
  {
    return 4;
  }
  int frictionConstraintSize_() const
  {
    return 2;
  }
  void initializeImpactAwareCopConstraint_();
  void computeImpactAwareCopConstraint_();
  void updateCopLogs_();

  void initializeImpactAwareFrictionConstraint_();
  void computeImpactAwareFrictionConstraint_();

  std::string mcContactName_;

  // Predictor
  mi_qpEstimator & predictor_;

  sva::ForceVecd measuredWrench_;

  // Alpha vector
  Eigen::VectorXd robotJointVelocity_;

  const ImpactAwareConstraintParams<Eigen::Vector2d> & constraintParams_;

  // dt * J_deltatau / impact_duration
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  Eigen::MatrixXd A_cop_;
  Eigen::VectorXd b_cop_;

  Eigen::MatrixXd A_friction_;
  Eigen::VectorXd b_friction_;


  Eigen::MatrixXd G_cop_;

  // Cop
  Eigen::Vector3d cop_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d cop_stateJump_ = Eigen::Vector3d::Zero();
  //Eigen::Vector2d cop_perturb_whole_ = Eigen::Vector2d::Zero();

  // Update sensor measurements. 
  


  // For the friction constraint
  Eigen::MatrixXd multiplier_; ///< (I - (1 + mu)nnT)
};

} // namespace mc_impact
