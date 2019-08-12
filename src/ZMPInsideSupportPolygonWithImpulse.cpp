#include "ZMPInsideSupportPolygonWithImpulse.h"

#include <mc_prediction/mi_impactPredictor.h>

namespace mc_impact
{

ZMPInsideSupportPolygonWithImpulse(const mc_solver::QPSolver & solver,
                                   const mc_rbdyn::Contact & leftSoleContact,
                                   const std::string & clName,
                                   const mc_rbdyn::Contact & rightSoleContact,
                                   const std::string & crName,
                                   mi_impactPredictor & predictor)
: InequalityConstraint(leftSoleContact.contactId(solver.robots()), rightSoleContact.contactId(solver.robots())),
  predictor_(predictor)
{

  // I am not sure how to fill in the selection matrix...

  slName_ = clName;
  srName_ = crName;
} // end of constructor

void ZMPInsideSupportPolygonWithImpulse::computeAb()
{

  // We need to calculate the support polygon:
  Eigen::MatrixXd A_s;
  A_s.resize(4, 6);
  A_s.setZero();

  A_s(0, 5) = -0.12;
  A_s(0, 1) = -1;
  A_s(1, 5) = -0.12;
  A_s(1, 1) = 1;
  A_s(2, 5) = -0.06;
  A_s(2, 0) = 1;
  A_s(3, 5) = -0.06;
  A_s(3, 0) = -1;

  A_ = A_s * (F_QP_left_sole_COM + F_QP_right_sole_COM);

  b_ = -A_s
       * (predictor_.getImpulsiveForceCOM("l_sole").vector() + predictor_.getImpulsiveForceCOM("r_sole").vector()
          + predictor_.getImpulsiveForceCOM().vector());

} // end of computeAb

} // namespace mc_impact
