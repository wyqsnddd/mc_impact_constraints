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

#include "BoundJointVelocityJump.h"

namespace mc_impact
{

BoundJointVelocityJump::BoundJointVelocityJump(mi_qpEstimator & predictor, double dt, double multi, bool debug)
: BoundJointVelocityJump(predictor,
                         dt,
                         multi * rbd::dofToVector(predictor.getSimRobot().mb(), predictor.getSimRobot().vl()),
                         multi * rbd::dofToVector(predictor.getSimRobot().mb(), predictor.getSimRobot().vu()),
                         debug)
{
}

BoundJointVelocityJump::BoundJointVelocityJump(mi_qpEstimator & predictor,
                                               double dt,
                                               const Eigen::VectorXd & LBound,
                                               const Eigen::VectorXd & UBound,
                                               bool debug)
: mc_solver::GenInequalityConstraint(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt),
  alpha_L_(LBound), alpha_U_(UBound), debug_(debug)
{
  alpha_.resize(alpha_L_.size());
  if(predictor_.getSimRobot().mb().joint(0).dof() == 6)
  {
    startIndex_ = 6;
  }
  int nDof = predictor_.getSimRobot().mb().nrDof();
  int realDof = nDof - startIndex_;
  alpha_L_ = alpha_L_.tail(realDof).eval();
  alpha_U_ = alpha_U_.tail(realDof).eval();
  A_.resize(realDof, nDof);
  L_.resize(realDof);
  U_.resize(realDof);

  diff_lower_.resize(realDof);
  diff_upper_.resize(realDof);
  test_delta_vel_.resize(realDof);
}

int BoundJointVelocityJump::maxGenInEq() const
{
  return static_cast<int>(alpha_L_.size());
}

void BoundJointVelocityJump::computeALU()
{
  const auto & robot = predictor_.getSimRobot();
  const auto & J_delta = predictor_.getJacobianDeltaAlpha();
  A_ = J_delta.block(startIndex_, 0, J_delta.rows() - startIndex_, J_delta.cols()) * dt_;
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  L_ = alpha_L_ - alpha_.tail(alpha_L_.size())
       - J_delta.block(startIndex_, 0, J_delta.rows() - startIndex_, J_delta.cols()) * alpha_;
  U_ = L_ - alpha_L_ + alpha_U_;

  if(debug_)
  {

    Eigen::VectorXd alphaD = (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD));

    diff_upper_ = U_ - A_ * alphaD;
    diff_lower_ = U_ - diff_upper_ - L_;

    test_delta_vel_ = J_delta * (alphaD * dt_ + alpha_);
  }
}

} // namespace mc_impact
