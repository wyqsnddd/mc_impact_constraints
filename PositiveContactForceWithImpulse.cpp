#include "PositiveContactForceWithImpulse.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

PositiveContactForceWithImpulse::PositiveContactForceWithImpulse(const mc_solver::QPSolver & solver,
                                                                 const mc_rbdyn::Contact & contact,
                                                                 mi_impactPredictor & predictor)
: InequalityConstraint(contact.contactId(solver.robots())),
  predictor_(predictor)
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
}

const Eigen::VectorXd & PositiveContactForceWithImpulse::bInEq() const
{
  return predictor_.getImpulsiveForce(sName_);
}


}
