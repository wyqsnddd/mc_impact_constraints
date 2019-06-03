#pragma once

#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ContactWrenchMatrixToLambdaMatrix.h"

/** Forward declaration **/
class mi_impactPredictor;

namespace mc_impact
{

/** This constraint restricts the ZMP inside the support polygon under impact.
 *
 * The impulsvie force at the impact body as well as the induced impulsive forces from the established contacts are
 * taken into account. The equivalent wrench of all of them at the COM \f$F_{ex}\f$ is used to predict the jump of the
 * ZMP under impact.
 *
 * Assuming the feet has planar contact surface, we use the convex hull enclosing the two foot contact areas to
 * construct the support polygon: \f[ A_{ZMP}F_{QP} \leqslant -A_{ZMP}F_{ex}. \f}
 */

struct ZMPInsideSupportPolygonWithImpulse : public mc_solver::InequalityConstraint
{
  ZMPInsideSupportPolygonWithImpulse(const mc_solver::QPSolver & solver,
                                     const mc_rbdyn::Contact & leftSoleContact,
                                     const std::string & clName,
                                     const mc_rbdyn::Contact & rightSoleContact,
                                     const std::string & crName,
                                     mi_impactPredictor & predictor);

  inline int maxInEq() const override
  {
    return 3;
  }

  inline std::string nameInEq() const override
  {
    return "ZMPInsideSupportPolygonWithImpulse";
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
  Eigen::MatrixXd multiplier_;
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;
  std::string slName_;
  std::string srName_;

}; // end of struct
} // namespace mc_impact
