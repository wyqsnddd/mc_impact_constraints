#include "COPInsideContactAreaWithImpulse.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

COPInsideContactAreaWithImpulse::COPInsideContactAreaWithImpulse(const mc_solver::QPSolver & solver,
                                                                 const mc_rbdyn::Contact & contact,
                                                                 mi_impactPredictor & predictor,
                                                                 double mu)
: InequalityConstraint(contact.contactId(solver.robots())),
  predictor_(predictor)
{
  auto cid = contact.contactId(solver.robots());
  ContactWrenchMatrixToLambdaMatrix transformer(solver, cid);
  Eigen::Vector3d normal = contact.X_0_r2s(solver.robots().robot(contact.r2Index())).rotation().row(2).transpose();
  multiplier_ = (Eigen::MatrixXd::Identity(3,3) - (1 + mu)*normal.transpose()*normal);
  Eigen::MatrixXd selector = Eigen::MatrixXd::Zero(3, 6);
  selector(0, 3) = 1;
  selector(1, 4) = 1;
  selector(2, 5) = 1;
  // selector * F = f => multiplier_ * selector * F = multiplier_ * f
  A_ = transformer.transform(multiplier_ * selector);
  sName_ = contact.r1Surface()->name();
}

void COPInsideContactAreaWithImpulse::computeAb()
{
  b_ = -multiplier_*predictor_.getImpulsiveForce(sName_);
}


}
