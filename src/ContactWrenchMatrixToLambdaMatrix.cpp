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

#include "ContactWrenchMatrixToLambdaMatrix.h"

#include <mc_rtc/logging.h>

namespace mc_impact
{

ContactWrenchMatrixToLambdaMatrix::ContactWrenchMatrixToLambdaMatrix(const mc_solver::QPSolver & solver,
                                                                     const tasks::qp::ContactId & id)
{
  auto qp_c = solver.contactById(id);
  if(qp_c.first == -1)
  {
    LOG_ERROR_AND_THROW(std::runtime_error,
                        "[ContactWrenchMatrixToLambdaMatrix] Provided contact has not been added to the solver yet")
  }
  const auto & contact = qp_c.second;
  const auto & cones = contact.r1Cones;
  const auto & points = contact.r1Points;
  size_t nLambda = 0;
  Eigen::MatrixXd Gi = Eigen::MatrixXd::Zero(6, 0);
  for(const auto & c : cones)
  {
    nLambda += c.generators.size();
    Gi.resize(3, std::max<int>(Gi.cols(), c.generators.size()));
  }
  transform_.resize(6, nLambda);
  Eigen::MatrixXd X_cf_pi_T = Eigen::MatrixXd::Zero(6, 3);
  int col = 0;
  for(size_t i = 0; i < points.size(); ++i)
  {
    sva::PTransformd X_b_pi{contact.r1Points[i]};
    sva::PTransformd X_cf_pi = X_b_pi * contact.X_b1_cf.inv();
    X_cf_pi_T.block(0, 0, 3, 3) = sva::vector3ToCrossMatrix(X_cf_pi.translation()) * X_cf_pi.rotation().transpose();
    X_cf_pi_T.block(3, 0, 3, 3) = X_cf_pi.rotation().transpose();
    for(size_t j = 0; j < cones[i].generators.size(); ++j)
    {
      Gi.col(j) = cones[i].generators[j];
    }
    transform_.block(0, col, 6, cones[i].generators.size()) = X_cf_pi_T * Gi.block(0, 0, 3, cones[i].generators.size());
    col += cones[i].generators.size();
  }
}

} // namespace mc_impact
