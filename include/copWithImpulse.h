#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

namespace mc_impact
{

struct newCoPArea
{
  double min_x;
  double max_x;
  double min_y;
  double max_y;
};

struct copWithImpulse : public mc_solver::InequalityConstraint
{
  copWithImpulse(mi_qpEstimator & predictor,
                 const std::string & bodyName,
                 const std::string & sensorName,
                 double dt,
                 double impact_dt,
                 // const mc_rbdyn::Contact & contact,
                 const newCoPArea & area);

  inline int maxInEq() const override
  {
    return 4;
  }

  inline std::string nameInEq() const override
  {
    return "copWithImpulse";
  }
  inline const Eigen::Vector2d & getCop() const 
  {
    return cop_;
  }
  inline const Eigen::Vector2d & getCopPerturb() const 
  {
    return cop_perturb_;
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

  void computeAb() override;

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
  Eigen::MatrixXd A_cop_;

  // Cop 
  Eigen::Vector2d cop_ = Eigen::Vector2d::Zero();
  Eigen::Vector2d cop_perturb_= Eigen::Vector2d::Zero();
  //Eigen::Vector2d cop_perturb_whole_ = Eigen::Vector2d::Zero();
};

} // namespace mc_impact
