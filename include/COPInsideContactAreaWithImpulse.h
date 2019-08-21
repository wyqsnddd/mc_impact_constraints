#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ContactWrenchMatrixToLambdaMatrix.h"

namespace mc_impact
{

struct CoPArea
{
  double min_x;
  double max_x;
  double min_y;
  double max_y;
};

/** This constraint enforces that the CoP within the contact area of an established contact,
 * i.e. it enforces \f$ A_{cop} (F_{QP} + F_{impact}) \leqslant 0 \f$ or \f$ A_{cop}F_{QP} \leqslant - A_{cop}F_{impact}
 * \f$.
 *
 * Note that \f$ F_{impact}\f$ denotes the sum of the impulsive reaction force and the equivalent wrenchs from the other
 * end-effectors.
 */
struct COPInsideContactAreaWithImpulse : public mc_solver::InequalityConstraint
{
  COPInsideContactAreaWithImpulse(const mc_solver::QPSolver & solver,
                                  const mc_rbdyn::Contact & contact,
                                  const CoPArea & area,
                                  mi_qpEstimator & predictor,
                                  const std::string & sName);

  inline int maxInEq() const override
  {
    return 4;
  }

  inline std::string nameInEq() const override
  {
    return "COPInsideContactAreaWithImpulse";
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
  Eigen::MatrixXd A_;
  Eigen::MatrixXd A_cop_;
  Eigen::VectorXd b_;
  std::string sName_;
};

} // namespace mc_impact
