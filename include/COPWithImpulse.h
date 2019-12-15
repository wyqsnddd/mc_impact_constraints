#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ConstraintUtils.h"

namespace mc_impact
{

struct COPWithImpulse : public mc_solver::InequalityConstraintRobot
{
  COPWithImpulse(const std::string & mcContactName,
                 mi_qpEstimator & predictor,
                 const ImpactAwareConstraintParams<Eigen::Vector2d> & impactAwareConstraintParams);

  inline int maxInEq() const override
  {
    return 4;
  }

  inline std::string nameInEq() const override
  {
    return "COPWithImpulse";
  }
  inline const Eigen::Vector2d & getCop() const
  {
    return cop_;
  }
  inline const Eigen::Vector2d & getCopPerturb() const
  {
    return cop_perturb_;
  }

  inline const Eigen::Vector2d & getCopPerturbWhole() const
  {
    return cop_perturb_whole_;
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

  // Alpha vector
  Eigen::VectorXd alpha_;

  const ImpactAwareConstraintParams<Eigen::Vector2d> & constraintParams_;
  // dt * J_deltatau / impact_duration
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;
  Eigen::MatrixXd A_cop_;

  // Cop
  Eigen::Vector2d cop_ = Eigen::Vector2d::Zero();
  Eigen::Vector2d cop_perturb_ = Eigen::Vector2d::Zero();
  Eigen::Vector2d cop_perturb_whole_ = Eigen::Vector2d::Zero();
};

} // namespace mc_impact
