#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/GenInequalityConstraint.h>

namespace mc_impact
{

/** This class imlements equation (25) of the Humanoids 2019 paper
 *
 * alpha_lower - alpha - J_dq * alpha <= J_dq * alphadD / dt <= alpha_upper - alpha - J_dq * alpha
 */
struct BoundJointVelocityJump : public mc_solver::GenInequalityConstraintRobot
{
  /** Use the robot module provided velocity bounds */
  BoundJointVelocityJump(mi_qpEstimator & predictor, double dt, double multiplier = 1.0, bool debu = false);

  BoundJointVelocityJump(mi_qpEstimator & predictor,
                         double dt,
                         const Eigen::VectorXd & LBound,
                         const Eigen::VectorXd & UBound,
                         bool debug);

  int maxGenInEq() const override;

  inline std::string nameGenInEq() const override
  {
    return "BoundJointVelocityJump";
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
    return diff_upper_;
  }
  inline const Eigen::VectorXd & getDiffLower()
  {
    return diff_lower_;
  }
  inline const Eigen::VectorXd & getDeltaVel()
  {
    return test_delta_vel_;
  }

  void compute() override;

private:
  // Predictor
  mi_qpEstimator & predictor_;
  // Timestep
  double dt_;
  // Lower joint velocity bound
  Eigen::VectorXd alpha_L_;
  // Upper joint velocity bound
  Eigen::VectorXd alpha_U_;
  // Alpha vector
  Eigen::VectorXd alpha_;
  // 0 for fixed-based, 6 otherwise
  int startIndex_ = 0;

  // J_deltaq / dt
  Eigen::MatrixXd A_;
  // alpha_min - alpha - J_deltaq * alpha
  Eigen::VectorXd L_;
  // alpha_max - alpha - J_deltaq * alpha
  Eigen::VectorXd U_;

  bool debug_;
  Eigen::VectorXd diff_lower_;
  Eigen::VectorXd diff_upper_;
  Eigen::VectorXd test_delta_vel_;
};

} // namespace mc_impact
