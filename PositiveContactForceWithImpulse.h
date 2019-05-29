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
 * i.e. it enforces \f$ f_{QP} + f_{impact} \geqslant 0 \f$ or \f$ -f_{QP} \leqslant f_{impact}
 *
 */
struct PositiveContactForceWithImpulse : public mc_solver::InequalityConstraint
{
  PositiveContactForceWithImpulse(const mc_solver::QPSolver & solver,
                                  const mc_rbdyn::Contact & contact,
                                  mi_impactPredictor & predictor);

  inline int maxInEq() const override { return 3; }

  inline std::string nameInEq() const override { return "PositiveContactForceWithImpulse"; }

  inline const Eigen::MatrixXd & A() const override { return A_; }

  const Eigen::VectorXd & bInEq() const override;

  inline void computeAb() override {}
private:
  mi_impactPredictor & predictor_;
  Eigen::MatrixXd A_;
  std::string sName_;
};

}
