#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ConstraintUtils.h"
#include <McDynamicStability/McComArea.h>
#include <McDynamicStability/McDCMArea.h>
#include <McDynamicStability/McZMPArea.h>

namespace mc_impact
{
struct ImpactAwareZMPConstraint : public mc_solver::InequalityConstraintRobot
{
  ///< ZMP defined with a set of points
  ImpactAwareZMPConstraint(mi_qpEstimator & predictor,
			   std::shared_ptr<McContactSet> contactSetPtr,
                           const ImpactAwareConstraintParams<Eigen::Vector2d> & params);

  /*!
   * \brief returns  the number of rows in realtime such that the QP can adjust the constraint size.
   */
  inline int nrInEq() const override
  {
    if(getParams().updateMcZMPArea)
    {
      return getMcZMPArea()->getMaxNumVertex();
    }
    else
    {
      return static_cast<int>(getParams().zmpAreaVertexSet.size());
    }
  }
  /*! \brief returns the upper bound of the amount of constraints for the QP to reserve memory accordingly.
   */
  inline int maxInEq() const override
  {
    if(getParams().updateMcZMPArea)
    {
      return getMcZMPArea()->getMaxNumVertex();
    }
    else
    {
      return static_cast<int>(getParams().zmpAreaVertexSet.size());
    }
  }

  inline const std::shared_ptr<const mc_impact::McZMPArea<Eigen::Vector2d>> getMcZMPArea() const
  {
    if(enabledMcZMPArea())
    {
      return mcZMPAreaPtr_;
    }
    else
    {
      LOG_ERROR_AND_THROW(std::runtime_error, "Asking for McZMPArea, which  is not initialized and updated.");
    }
  }
  inline std::string nameInEq() const override
  {
    return "ImpactAwareZMPConstraint";
  }
  inline const Eigen::MatrixXd & A() const override
  {
    return A_;
  }
  inline const Eigen::VectorXd & bInEq() const override
  {
    return b_;
  }

  inline const ImpactAwareConstraintParams<Eigen::Vector2d> & getParams() const
  {
    return params_;
  }

  inline const FloatingBaseStates & getFloatingBaseStates() const
  {
    return floatingBaseStates_;
  }

  bool enabledMcZMPArea() const
  {
    return getParams().updateMcZMPArea;
  }

  void compute() override;

  inline const IeqConstraintBlocks & getIeqBlocksZMP() const
  {
    return ieqConstraintBlocksZMP_;
  }

  inline void setIeqBlocksZMP(const IeqConstraintBlocks & input)
  {
    ieqConstraintBlocksZMP_ = input;
  }
  /*! \return reference to the robot (either pyhsices-engine simulated or real)
   */

  inline const mc_rbdyn::Robot & robot() const
  {
    return robot_;
  }

  /*! \return if the 'samplePoint' is inside the multi-contact ZMP area.
   */
  bool pointInsideMcZMPArea(const Eigen::Vector3d & samplePoint) const;

  void printZMPConstraintMatrix() const;

  /*! \return The constant: \omega = sqrt(9.8/com-z)
   */
  inline double getOmega() const
  {
    return omega_;
  }

private:
  // Predictor
  mi_qpEstimator & predictor_;
  ImpactAwareConstraintParams<Eigen::Vector2d> params_;
  std::shared_ptr<McContactSet> contactSetPtr_;

  const mc_rbdyn::Robot & robot_;

  inline void calcOmega_(const double & c_z)
  {
    assert(c_z > 0.0);
    omega_ = sqrt(9.8 / c_z);
  }

  double omega_;


  inline int dof_() const
  {
    return robot_.mb().nrDof();
  }

  // Multi-contact ZMP area calculator:
  std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr_;

  void updateZMPConstraint_();

  void updateMcAreas_(double height);
  void updateMcZMPAreas_(double height);
  void getZMPBlocks(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench);
  void calculateZMP_();

  std::shared_ptr<rbd::CoMJacobian> comJacobianPtr_;

  // Left hand side of the ZMP and DCM constraint
  Eigen::MatrixXd A_zmp_;
  Eigen::MatrixXd A_dcm_;
  void fixedSupportPolygonSetup_();

  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  // Needs to be reset in each iteration if we update the Multi-contact areas in realtime.
  void updateAbMatrixsize_();

  IeqConstraintBlocks ieqConstraintBlocksZMP_;

    FloatingBaseStates floatingBaseStates_;
  void updateFloatingBaseState_();

  Eigen::VectorXd robotJointVelocity_;
}; // End of struct name: ImpactAwareZMPConstraint

} // end of namespace mc_impact
