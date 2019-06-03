#include "ZeroSlippageWithImpulse.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

ZeroSlippageWithImpulse::ZeroSlippageWithImpulse(const mc_solver::QPSolver & solver,
                                                                 const mc_rbdyn::Contact & contact,
                                                                 mi_impactPredictor & predictor,
                                                                 const std::string & cBody,
                                                                 double mu)
: InequalityConstraint(contact.contactId(solver.robots())), predictor_(predictor)
{
  auto cid = contact.contactId(solver.robots());
  ContactWrenchMatrixToLambdaMatrix transformer(solver, cid);
  Eigen::Vector3d normal = contact.X_0_r2s(solver.robots()).rotation().row(2).transpose();

  // Should we use normal*normal.transpose()?
  multiplier_ = (Eigen::MatrixXd::Identity(3, 3) - (1 + mu) * normal * normal.transpose());
  Eigen::MatrixXd selector = Eigen::MatrixXd::Zero(3, 6);
  selector(0, 3) = 1;
  selector(1, 4) = 1;
  selector(2, 5) = 1;
  // selector * F = f => multiplier_ * selector * F = multiplier_ * f
  A_ = transformer.transform(multiplier_ * selector);
  sName_ = cBody;
}

void ZeroSlippageWithImpulse::computeAb()
{
  const auto & force = predictor_.getImpulsiveForce(sName_);
  if(force.rows() == 3 && force.cols() == 1)
  {
    b_ = -multiplier_ * predictor_.getImpulsiveForce(sName_);
  }
  else
  {
    b_.setZero(3);
  }
}

} // namespace mc_impact
