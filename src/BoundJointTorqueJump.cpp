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

#include "BoundJointTorqueJump.h"

namespace mc_impact
{

BoundJointTorqueJump::BoundJointTorqueJump(mi_qpEstimator & predictor,
                                           double dt,
                                           double impact_dt,
                                           double mult,
                                           bool debug)
: BoundJointTorqueJump(predictor,
                       dt,
                       impact_dt,
                       mult * rbd::dofToVector(predictor.getSimRobot().mb(), predictor.getSimRobot().tl()),
                       mult * rbd::dofToVector(predictor.getSimRobot().mb(), predictor.getSimRobot().tu()),
                       debug)
{
}

BoundJointTorqueJump::BoundJointTorqueJump(mi_qpEstimator & predictor,
                                           double dt,
                                           double impact_dt,
                                           const Eigen::VectorXd & LBound,
                                           const Eigen::VectorXd & UBound,
                                           bool debug)
: mc_solver::GenInequalityConstraint(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt),
  impact_dt_(impact_dt), tau_L_(LBound), tau_U_(UBound), debug_(debug)
{
  // std::cout<<"The lower impulsive torque bound is: "<<std::endl<<LBound.transpose()<<std::endl;
  // std::cout<<"The upper impulsive torque bound is: "<<std::endl<<UBound.transpose()<<std::endl;
  alpha_.resize(tau_L_.size());
  if(predictor_.getSimRobot().mb().joint(0).dof() == 6)
  {
    startIndex_ = 6;
  }
  int nDof = predictor_.getSimRobot().mb().nrDof();
  int realDof = nDof - startIndex_;
  tau_L_ = tau_L_.tail(realDof).eval();
  tau_U_ = tau_U_.tail(realDof).eval();
  A_.resize(realDof, nDof);
  L_.resize(realDof);
  U_.resize(realDof);
  test_delta_torque_.resize(realDof);
  difference_lower_.resize(realDof);
  difference_upper_.resize(realDof);
}

int BoundJointTorqueJump::maxGenInEq() const
{
  return static_cast<int>(tau_L_.size());
}

void BoundJointTorqueJump::computeALU()
{
  const auto & robot = predictor_.getSimRobot();
  const auto & J_deltatau = predictor_.getJacobianDeltaTau();
  A_ = (dt_ / impact_dt_) * J_deltatau.block(startIndex_, 0, J_deltatau.rows() - startIndex_, J_deltatau.cols());
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  L_ = tau_L_
       - J_deltatau.block(startIndex_, 0, J_deltatau.rows() - startIndex_, J_deltatau.cols()) * alpha_ / impact_dt_;
  U_ = L_ - tau_L_ + tau_U_;

  if(debug_)
  {

    Eigen::VectorXd alphaD = (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD));

    difference_upper_ = U_ - A_ * alphaD;
    difference_lower_ = U_ - difference_upper_ - L_;

    test_delta_torque_ = (1.0 / impact_dt_) * J_deltatau * (alpha_ + alphaD * dt_);
  }
}

} // namespace mc_impact
