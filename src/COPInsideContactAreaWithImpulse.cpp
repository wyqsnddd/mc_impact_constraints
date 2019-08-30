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

#include "COPInsideContactAreaWithImpulse.h"

namespace mc_impact
{

COPInsideContactAreaWithImpulse::COPInsideContactAreaWithImpulse(const mc_solver::QPSolver & solver,
                                                                 const mc_rbdyn::Contact & contact,
                                                                 const CoPArea & area,
                                                                 mi_qpEstimator & predictor,
                                                                 const std::string & sName)
: InequalityConstraint(contact.contactId(solver.robots())), predictor_(predictor)
{
  auto cid = contact.contactId(solver.robots());
  ContactWrenchMatrixToLambdaMatrix transformer(solver, cid);
  sName_ = sName;

  A_cop_ = Eigen::MatrixXd::Zero(4, 6);

  // - cy - max_x * fz
  A_cop_(0, 1) = -1;
  A_cop_(0, 5) = -area.max_x;

  // cy + min_x * fz
  A_cop_(1, 1) = 1;
  A_cop_(1, 5) = area.min_x;

  // cx - max_y * fz
  A_cop_(2, 0) = 1;
  A_cop_(2, 5) = -area.max_y;

  // - cx + min_y * fz
  A_cop_(3, 0) = -1;
  A_cop_(3, 5) = area.min_y;

  A_ = transformer.transform(A_cop_);
}

void COPInsideContactAreaWithImpulse::computeAb()
{
  Eigen::Vector6d f;
  f.setZero();
  f.segment(3, 3) = predictor_.getEndeffector(sName_).estimatedAverageImpulsiveForce;

  b_.noalias() = -A_cop_ * f;
}

} // namespace mc_impact
