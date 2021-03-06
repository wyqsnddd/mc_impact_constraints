/* Copyright 2018-2019 CNRS-UM LIRMM
 *
 * \author Yuquan Wang, Arnaud Tanguy and Pierre Gergondet
 *
 * 
 *
 * mc_impact_constraints is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * mc_impact_constraints is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with mc_impact_constraints. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

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
                 bool allforce = true,
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
  inline const Eigen::Vector3d & getZMP_prediction_feetforce()
  {
    if(debug_)
      return zmpPrediction_feetforce_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_prediction_allforce()
  {
    if(debug_)
      return zmpPrediction_allforce_;
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
  const ZMPArea area_;
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
  bool allForce_;

  void calcZMP_();
  Eigen::Vector3d zmpSensor_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPerturbation_ = Eigen::Vector3d::Zero();

  Eigen::Vector3d zmpPrediction_allforce_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPrediction_feetforce_ = Eigen::Vector3d::Zero();

  Eigen::Vector4d difference_ = Eigen::Vector4d::Zero();
};

} // namespace mc_impact
