#include "COPWithImpulse.h"

namespace mc_impact
{

COPWithImpulse::COPWithImpulse(mi_qpEstimator & predictor,
                               const std::string & bodyName,
                               const std::string & sensorName,
                               double dt,
                               double impact_dt,
                               // const mc_rbdyn::Contact & contact,
                               const newCoPArea & area)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt), impact_dt_(impact_dt),
  area_(area)
{

  bName_ = bodyName;
  sName_ = sensorName;

  A_cop_ = Eigen::MatrixXd::Zero(4, 6);

  // - cy - max_x * fz
  A_cop_(0, 1) = -1;
  A_cop_(0, 5) = -area_.max_x;

  // cy + min_x * fz
  A_cop_(1, 1) = 1;
  A_cop_(1, 5) = area_.min_x;

  // cx - max_y * fz
  A_cop_(2, 0) = 1;
  A_cop_(2, 5) = -area_.max_y;

  // - cx + min_y * fz
  A_cop_(3, 0) = -1;
  A_cop_(3, 5) = area_.min_y;

  int nDof = predictor_.getSimRobot().mb().nrDof();

  alpha_.resize(nDof);
  A_.resize(4, nDof);
  b_.resize(4);
}

void COPWithImpulse::compute()
{

  const auto & robot = predictor_.getSimRobot();
  const auto & J_deltaF = predictor_.getJacobianDeltaF(bName_);
  // std::cout<<"size of J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;
  // std::cout<<"size of reduced J_deltaF: "<<J_deltaF.rows()<<", "<<J_deltaF.cols()<<std::endl;

  // A_ = (dt_ / impact_dt_) * A_cop_.block(0, 3, 4, 3)*J_deltaF.block(0, startIndex_, J_deltaF.rows(), J_deltaF.cols()
  // - startIndex_);
  A_ = (dt_ / impact_dt_) * A_cop_.block(0, 3, 4, 3) * J_deltaF;
  // std::cout<<"size of A_: "<<A_.rows()<<", "<<A_.cols()<<std::endl;
  rbd::paramToVector(robot.mbc().alpha, alpha_);

  // std::cout<<"size of alpha_"<<alpha_.rows()<<std::endl;
  sva::ForceVecd bodyWrenchSensor = predictor_.getSimRobot().bodyWrench(bName_);
  b_ = -(A_cop_ * bodyWrenchSensor.vector() + A_cop_.block(0, 3, 4, 3) * J_deltaF * alpha_ / impact_dt_);
  /*
  b_ = -(A_cop_ * predictor_.getSimRobot().forceSensor(sName_).wrench().vector()
         + A_cop_.block(0, 3, 4, 3) * J_deltaF * alpha_ / impact_dt_);
   */
  cop_.x() = -bodyWrenchSensor.couple().y() / bodyWrenchSensor.force().z();
  cop_.y() = bodyWrenchSensor.couple().x() / bodyWrenchSensor.force().z();

  double perturbedNormalForce =
      bodyWrenchSensor.force().z() + predictor_.getEndeffector(bName_).estimatedAverageImpulsiveForce.z();

  cop_perturb_.x() = -bodyWrenchSensor.couple().y() / perturbedNormalForce;
  cop_perturb_.y() = bodyWrenchSensor.couple().x() / perturbedNormalForce;

  cop_perturb_whole_.setZero();
  // std::cout<<"copConstraint: The converted wrench of "<<bName_<<" is: "<<
  // predictor_.getEndeffector(bName_).perturbedWrench.vector().transpose()<<std::endl;
  double perturbedNormalForce_whole =
      perturbedNormalForce + predictor_.getEndeffector(bName_).perturbedWrench.force().z();

  cop_perturb_whole_.x() =
      -(bodyWrenchSensor.couple().y() + predictor_.getEndeffector(bName_).perturbedWrench.couple().y())
      / perturbedNormalForce_whole;

  cop_perturb_whole_.y() =
      (bodyWrenchSensor.couple().x() + predictor_.getEndeffector(bName_).perturbedWrench.couple().x())
      / perturbedNormalForce_whole;
}

} // namespace mc_impact
