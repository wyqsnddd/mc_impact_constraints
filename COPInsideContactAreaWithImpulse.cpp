#include "COPInsideContactAreaWithImpulse.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

COPInsideContactAreaWithImpulse::COPInsideContactAreaWithImpulse(const mc_solver::QPSolver & solver,
                                                                 const mc_rbdyn::Contact & contact,
                                                                 const CoPArea & area,
                                                                 mi_impactPredictor & predictor,
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
  // Eigen::VectorXd f = predictor_.getImpulsiveForceCOP(sName_).vector();
  Eigen::Vector6d f;
  f.setZero();
  f.segment(3, 3) = predictor_.getImpulsiveForce(sName_);

  b_.noalias() = -A_cop_ * f;
}

} // namespace mc_impact
