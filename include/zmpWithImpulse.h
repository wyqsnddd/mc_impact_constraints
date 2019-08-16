#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

/** Forward declaration */
class mi_impactPredictor;

namespace mc_impact
{

struct ZMPArea
{
  double min_x;
  double max_x;
  double min_y;
  double max_y;
};

struct supportContact
{
  std::string bodyName;
  std::string sensorName;
};

struct zmpWithImpulse : public mc_solver::InequalityConstraint
{
  zmpWithImpulse(mi_qpEstimator & predictor,
                 const std::vector<supportContact> & supports,
                 double dt,
                 double impact_dt,
                 const ZMPArea & area,
                 bool debug = false);

  inline int maxInEq() const override
  {
    return 4;
  }

  inline std::string nameInEq() const override
  {
    return "zmpWithImpulse";
  }

  inline const Eigen::MatrixXd & A() const override
  {
    return A_;
  }
  inline const Eigen::VectorXd & bInEq() const override
  {
    return b_;
  }

  void computeAb() override;

  void getInertialItems(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench);
  inline const Eigen::MatrixXd & getA()
  {
    return A_;
  }
  inline const Eigen::VectorXd & getb()
  {
    return b_;
  }
  inline const Eigen::MatrixXd & getZMP()
  {
    return A_zmp_;
  }
  inline const Eigen::Vector3d & getZMP_sensor()
  {
    if(debug_)
      return zmpSensor_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_perturbation()
  {
    if(debug_)
      return zmpPerturbation_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_prediction()
  {
    if(debug_)
      return zmpPrediction_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector4d & getZMP_constraint_difference()
  {
    if(debug_)
      return difference_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }

private:
  // Predictor
  mi_qpEstimator & predictor_;
  // Timestep
  double dt_;
  // Impact duration
  double impact_dt_;
  // Alpha vector
  Eigen::VectorXd alpha_;
  std::vector<supportContact> supports_;
  // Name of the body of interest.
  // std::string bName_;
  // Name of the force sensor
  // std::string sName_;
  // dt * J_deltatau / impact_duration
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;
  Eigen::MatrixXd A_zmp_;
  Eigen::MatrixXd A_zmp_f_;
  bool debug_;
  void calcZMP_(const Eigen::MatrixXd & sumJac, const Eigen::Vector6d & exWrench);
  Eigen::Vector3d zmpSensor_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPerturbation_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPrediction_ = Eigen::Vector3d::Zero();
  Eigen::Vector4d difference_ = Eigen::Vector4d::Zero();
};

} // namespace mc_impact
