#pragma once

#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ContactWrenchMatrixToLambdaMatrix.h"

/** Forward declaration */
class mi_impactPredictor;

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
                                  mi_impactPredictor & predictor,
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
  mi_impactPredictor & predictor_;
  Eigen::MatrixXd multiplier_; // (I - (1 + mu)nnT)
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;
  std::string sName_;
};

} // namespace mc_impact
