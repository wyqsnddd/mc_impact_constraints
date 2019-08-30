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

#include "ZeroSlippageWithImpulse.h"

namespace mc_impact
{

ZeroSlippageWithImpulse::ZeroSlippageWithImpulse(const mc_solver::QPSolver & solver,
                                                 const mc_rbdyn::Contact & contact,
                                                 mi_qpEstimator & predictor,
                                                 const std::string & cBody,
                                                 double mu)
: InequalityConstraint(contact.contactId(solver.robots())), predictor_(predictor)
{
  auto cid = contact.contactId(solver.robots());
  ContactWrenchMatrixToLambdaMatrix transformer(solver, cid);
  Eigen::Vector3d normal = contact.X_0_r2s(solver.robots()).rotation().row(2).transpose();

  multiplier_ = (Eigen::MatrixXd::Identity(3, 3) - (1 + mu) * normal * normal.transpose());
  Eigen::MatrixXd selector = Eigen::MatrixXd::Zero(3, 6);
  selector(0, 3) = 1;
  selector(1, 4) = 1;
  selector(2, 5) = 1;
  A_ = transformer.transform(multiplier_ * selector);
  sName_ = cBody;
}

void ZeroSlippageWithImpulse::computeAb()
{
  const auto & force = predictor_.getEndeffector(sName_).estimatedAverageImpulsiveForce;
  if(force.rows() == 3 && force.cols() == 1)
  {
    b_ = -multiplier_ * predictor_.getEndeffector(sName_).estimatedAverageImpulsiveForce;
  }
  else
  {
    b_.setZero(3);
  }
}

} // namespace mc_impact
