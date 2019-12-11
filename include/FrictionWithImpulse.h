#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

namespace mc_impact
{

/** This class imlements equation (25) of the Humanoids 2019 paper
 *
 */
struct FrictionWithImpulse : public mc_solver::InequalityConstraintRobot
{
  FrictionWithImpulse(mi_qpEstimator & predictor,
                      const std::string & bodyName,
                      const std::string & sensorName,
                      //const mc_rbdyn::Contact & contact,
                      double dt,
                      double impact_dt,
                      double mu = mc_rbdyn::Contact::defaultFriction);

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

private:
  // Predictor
  mi_qpEstimator & predictor_;
  // Timestep
  double dt_;
  // Impact duration
  double impact_dt_;
  // Alpha vector
  Eigen::VectorXd alpha_;
  // Name of the body of interest.
  std::string bName_;
  // Name of the force sensor
  std::string sName_;
  // dt * J_deltatau / impact_duration
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  Eigen::MatrixXd multiplier_; // (I - (1 + mu)nnT)
};

} // namespace mc_impact
