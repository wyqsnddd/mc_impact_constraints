#include "COPWithImpulse.h"

namespace mc_impact
{

COPWithImpulse::COPWithImpulse(mi_qpEstimator & predictor,
                               double dt,
                               double impactDuration,
                               const McContactParams & contactParams)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt),
  impactDuration_(impactDuration), mcContactParams_(contactParams)
{

  A_cop_ = Eigen::MatrixXd::Zero(4, 6);

  double maxX = getParams().halfX;
  double minX = -getParams().halfX;
  double maxY = getParams().halfY;
  double minY = -getParams().halfY;

  // - cy - max_x * fz
  A_cop_(0, 1) = -1;
  A_cop_(0, 5) = -maxX;

  // cy + min_x * fz
  A_cop_(1, 1) = 1;
  A_cop_(1, 5) = minX;

  // cx - max_y * fz
  A_cop_(2, 0) = 1;
  A_cop_(2, 5) = -maxY;

  // - cx + min_y * fz
  A_cop_(3, 0) = -1;
  A_cop_(3, 5) = minY;

  int nDof = predictor_.getSimRobot().mb().nrDof();

  alpha_.resize(nDof);
  A_.resize(4, nDof);
  b_.resize(4);
}

void COPWithImpulse::compute()
{

  const auto & robot = predictor_.getSimRobot();
  const auto & J_deltaF = predictor_.getJacobianDeltaF(getParams().bodyName);
  // std::cout<<"size of J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;
  // std::cout<<"size of reduced J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;

  // A_ = (dt_ / impact_dt_) * A_cop_.block(0, 3, 4, 3)*J_deltaF.block(0, startIndex_, J_deltaF.rows(), J_deltaF.cols()
  // - startIndex_);
  A_ = (dt_ / impactDuration_) * A_cop_.block(0, 3, 4, 3) * J_deltaF;
  // std::cout<<"size of A_: "<<A_.rows()<<", "<<A_.cols()<<std::endl;
  rbd::paramToVector(robot.mbc().alpha, alpha_);

  // std::cout<<"size of alpha_"<<alpha_.rows()<<std::endl;
  sva::ForceVecd bodyWrenchSensor = predictor_.getSimRobot().bodyWrench(getParams().bodyName);
  b_ = -(A_cop_ * bodyWrenchSensor.vector() + A_cop_.block(0, 3, 4, 3) * J_deltaF * alpha_ / impactDuration_);
  /*
  b_ = -(A_cop_ * predictor_.getSimRobot().forceSensor(sName_).wrench().vector()
         + A_cop_.block(0, 3, 4, 3) * J_deltaF * alpha_ / impact_dt_);
   */
  cop_.x() = -bodyWrenchSensor.couple().y() / bodyWrenchSensor.force().z();
  cop_.y() = bodyWrenchSensor.couple().x() / bodyWrenchSensor.force().z();

  double perturbedNormalForce =
      bodyWrenchSensor.force().z() + predictor_.getEndeffector(getParams().bodyName).estimatedAverageImpulsiveForce.z();

  cop_perturb_.x() = -bodyWrenchSensor.couple().y() / perturbedNormalForce;
  cop_perturb_.y() = bodyWrenchSensor.couple().x() / perturbedNormalForce;

  cop_perturb_whole_.setZero();
  // std::cout<<"copConstraint: The converted wrench of "<<bName_<<" is: "<<
  // predictor_.getEndeffector(bName_).perturbedWrench.vector().transpose()<<std::endl;
  double perturbedNormalForce_whole =
      perturbedNormalForce + predictor_.getEndeffector(getParams().bodyName).perturbedWrench.force().z();

  cop_perturb_whole_.x() =
      -(bodyWrenchSensor.couple().y() + predictor_.getEndeffector(getParams().bodyName).perturbedWrench.couple().y())
      / perturbedNormalForce_whole;

  cop_perturb_whole_.y() =
      (bodyWrenchSensor.couple().x() + predictor_.getEndeffector(getParams().bodyName).perturbedWrench.couple().x())
      / perturbedNormalForce_whole;
}

} // namespace mc_impact
