#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/GenInequalityConstraint.h>

namespace mc_impact
{

/** This class imlements equation (25) of the Humanoids 2019 paper
 *
 * min_deltau - J_deltatau * alpha / impact_duration <= dt * J_deltatau * alphaD / impact_duration
 *
 * dt * J_deltatau * alphaD / impact_duration <= max_deltau - J_deltatau * alpha / impact_duration
 */
struct BoundJointTorqueJump : public mc_solver::GenInequalityConstraint
{
  /** Multiply the regular torque bounds limits by the provided multiplier */
  BoundJointTorqueJump(mi_qpEstimator & predictor, double dt, double impact_dt, double mult, bool debug = false);

  BoundJointTorqueJump(mi_qpEstimator & predictor,
                       double dt,
                       double impact_dt,
                       const Eigen::VectorXd & LBound,
                       const Eigen::VectorXd & UBound,
                       bool debug = false);

  int maxGenInEq() const override;

  inline std::string nameGenInEq() const override
  {
    return "BoundJointTorqueJump";
  }

  inline const Eigen::VectorXd & LowerGenInEq() const override
  {
    return L_;
  }
  inline const Eigen::VectorXd & UpperGenInEq() const override
  {
    return U_;
  }
  inline const Eigen::MatrixXd & A() const override
  {
    return A_;
  }
  inline const Eigen::VectorXd & getDiffUpper()
  {
    return difference_upper_;
  }
  inline const Eigen::VectorXd & getDiffLower()
  {
    return difference_lower_;
  }

  inline const Eigen::VectorXd & getDeltaTau()
  {
    return test_delta_torque_;
  }
private:
  void computeALU() override;

  // Predictor
  mi_qpEstimator & predictor_;
  // Timestep
  double dt_;
  // Impact duration
  double impact_dt_;
  // Lower impact torque bounds
  Eigen::VectorXd tau_L_;
  // Upper impact torque bounds
  Eigen::VectorXd tau_U_;
  // Alpha vector
  Eigen::VectorXd alpha_;
  // 0 for fixed-based, 6 otherwise
  int startIndex_ = 0;

  bool debug_;
  Eigen::VectorXd difference_upper_;
  Eigen::VectorXd difference_lower_;
  Eigen::VectorXd test_delta_torque_;

  // dt * J_deltatau / impact_duration
  Eigen::MatrixXd A_;
  // min_deltatau - J_deltatau * alpha / impact_duration
  Eigen::VectorXd L_;
  // max_deltatau - J_deltatau * alpha / impact_duration
  Eigen::VectorXd U_;
};

} // namespace mc_impact
