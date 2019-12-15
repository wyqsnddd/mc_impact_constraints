#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ConstraintUtils.h"

namespace mc_impact
{

/** This class imlements equation (25) of the Humanoids 2019 paper
 *
 */
struct FrictionWithImpulse : public mc_solver::InequalityConstraintRobot
{
  FrictionWithImpulse(const std::string & mcContactName, mi_qpEstimator & predictor, const ImpactAwareConstraintParams<Eigen::Vector2d> & impactAwareConstraintParams);

  inline int maxInEq() const override
  {
    return 2;
  }

  inline std::string nameInEq() const override
  {
    return "FrictionWithImpulse";
  }

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

private:

  std::string mcContactName_;

  // Predictor
  mi_qpEstimator & predictor_;
  
  const ImpactAwareConstraintParams<Eigen::Vector2d> & constraintParams_;
   // Alpha vector
  Eigen::VectorXd alpha_;
  // dt * J_deltatau / impact_duration
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  Eigen::MatrixXd multiplier_; // (I - (1 + mu)nnT)
};

} // namespace mc_impact
