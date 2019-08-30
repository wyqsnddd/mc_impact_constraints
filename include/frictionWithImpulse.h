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

/** This class imlements equation (25) of the Humanoids 2019 paper
 *
 */
struct frictionWithImpulse : public mc_solver::InequalityConstraint
{
  frictionWithImpulse(mi_qpEstimator & predictor,
                      const std::string & bodyName,
                      const std::string & sensorName,
                      const mc_rbdyn::Contact & contact,
                      double dt,
                      double impact_dt,
                      double mu = mc_rbdyn::Contact::defaultFriction);

  inline int maxInEq() const override
  {
    return 2;
  }

  inline std::string nameInEq() const override
  {
    return "frictionWithImpulse";
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
