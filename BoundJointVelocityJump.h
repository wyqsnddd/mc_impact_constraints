#pragma once

#include <mc_solver/GenInequalityConstraint.h>

/** Forward declaration */
class mi_impactPredictor;

namespace mc_impact
{

/** This class imlements equation (25) of the Humanoids 2019 paper
 *
 * alpha_lower - alpha - J_dq * alpha <= J_dq * alphadD / dt <= alpha_upper - alpha - J_dq * alpha
 */
struct BoundJointVelocityJump : public mc_solver::GenInequalityConstraint
{
  /** Use the robot module provided velocity bounds */
  BoundJointVelocityJump(mi_impactPredictor & predictor, double dt);

  BoundJointVelocityJump(mi_impactPredictor & predictor, double dt, const Eigen::VectorXd & LBound, const Eigen::VectorXd & UBound);

  int maxGenInEq() const override;

  inline std::string nameGenInEq() const override { return "BoundJointVelocityJump"; }

  inline const Eigen::VectorXd & LowerGenInEq() const override { return L_; }
  inline const Eigen::VectorXd & UpperGenInEq() const override { return U_; }
  inline const Eigen::MatrixXd & A() const override { return A_; }
private:
  void computeALU() override;

  // Predictor
  mi_impactPredictor & predictor_;
  // Timestep
  double dt_;
  // Lower joint velocity bound
  Eigen::VectorXd alpha_L_;
  // Upper joint velocity bound
  Eigen::VectorXd alpha_U_;
  // Alpha vector
  Eigen::VectorXd alpha_;

  // J_deltaq / dt
  Eigen::MatrixXd A_;
  // alpha_min - alpha - J_deltaq * alpha
  Eigen::VectorXd L_;
  // alpha_max - alpha - J_deltaq * alpha
  Eigen::VectorXd U_;
};

}
