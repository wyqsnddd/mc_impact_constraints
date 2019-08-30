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

#include "ContactWrenchMatrixToLambdaMatrix.h"

namespace mc_impact
{

/** This constraint enforces that the sum of the QP force and the impulse remains positive
 *
 * For a given contact with the environment where the normal is given by \f$n\f$, this enforces:
 *
 * \f[
 *   (I - (1 + \mu)nn^{T})(f_{QP} + f_{impulse}) \leqslant 0
 * \f]
 *
 */
struct ZeroSlippageWithImpulse : public mc_solver::InequalityConstraint
{
  ZeroSlippageWithImpulse(const mc_solver::QPSolver & solver,
                          const mc_rbdyn::Contact & contact,
                          mi_qpEstimator & predictor,
                          const std::string & cBody,
                          double mu = mc_rbdyn::Contact::defaultFriction);

  inline int maxInEq() const override
  {
    return 3;
  }

  inline std::string nameInEq() const override
  {
    return "ZeroSlippageWithImpulse";
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
  mi_qpEstimator & predictor_;
  Eigen::MatrixXd multiplier_; // (I - (1 + mu)nnT)
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;
  std::string sName_;
};

} // namespace mc_impact
