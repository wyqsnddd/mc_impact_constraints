#include "COPInsideContactAreaWithImpulse.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

COPInsideContactAreaWithImpulse::COPInsideContactAreaWithImpulse(const mc_solver::QPSolver & solver,
                                                                 const mc_rbdyn::Contact & contact,
                                                                 mi_impactPredictor & predictor)
: InequalityConstraint(contact.contactId(solver.robots())), predictor_(predictor)
{
  auto cid = contact.contactId(solver.robots());
  ContactWrenchMatrixToLambdaMatrix transformer(solver, cid);
  /** \f$ selector * F = -f \f$ */
  Eigen::MatrixXd selector = Eigen::MatrixXd::Zero(3, 6);
  selector(0, 3) = -1;
  selector(1, 4) = -1;
  selector(2, 5) = -1;
  A_ = transformer.transform(selector);
  sName_ = contact.r1Surface()->name();

  /// Suppose the contact area is x \f$ x \in [-0.12, 0.12], y \in [-0.06, 0.06]\f$

  A_cop.resize(4, 6);
  A_cop.setZero();

  A_cop(0, 5) = -0.12;
  A_cop(0, 1) = -1;
  A_cop(1, 5) = -0.12;
  A_cop(1, 1) = 1;
  A_cop(2, 5) = -0.06;
  A_cop(2, 0) = 1;
  A_cop(3, 5) = -0.06;
  A_cop(3, 0) = -1;

}

const Eigen::VectorXd & COPInsideContactAreaWithImpulse::bInEq() const
{
  return A_cop* predictor_.getImpulsiveForceCOP(sName_).vector();
}

} // namespace mc_impact
