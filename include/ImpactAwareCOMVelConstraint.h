# pragma once 

#include <mc_prediction/mi_qpEstimator.h> 

#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ConstraintUtils.h"

namespace mc_impact
{
struct ImpactAwareCOMVelConstraint : public mc_solver::InequalityConstraintRobot
{

  ImpactAwareCOMVelConstraint(std::shared_ptr<mi_qpEstimator> predictorPtr,
		  const ImpactAwareConstraintParams<Eigen::Vector2d> & params);

  /*!
   * \brief returns  the number of rows in realtime such that the QP can adjust the constraint size.
   */
  inline int nrInEq() const override
  {
    // Should be four, x-y
    return 2; 
  }

  /*! \brief returns the upper bound of the amount of constraints for the QP to reserve memory accordingly.
   */
  inline int maxInEq() const override
  {
    return nrInEq(); 
  }
  inline std::string nameInEq() const override
  {
    return "ImpactAwareCOMVelConstraint";
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
  inline const COMStates & getCOMStates() const
  {
    return COMStates_;
  }

  void compute() override;

  inline const mc_rbdyn::Robot & robot() const
  {
    return robot_;
  }

  /*! \return The constant: \omega = sqrt(9.8/com-z)
   */
  inline double getOmega() const
  {
    return omega_;
  }

  inline std::shared_ptr<mi_qpEstimator> getPredictor()
  {
    return predictorPtr_; 
  }
  
protected:

  std::shared_ptr<mi_qpEstimator> predictorPtr_;

  ImpactAwareConstraintParams<Eigen::Vector2d> params_;
  std::shared_ptr<McContactSet> contactSetPtr_;

  const mc_rbdyn::Robot & robot_;

  inline void updateGains_()
  {
    calcOmega_(COMStates_.Com.z());

    pGain_ = 1 + poleOne_*poleTwo_; 
    dGain_ = -(poleOne_ + poleTwo_)/getOmega();
  }

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

  void updateCOMVelConstraint_();
  std::shared_ptr<rbd::CoMJacobian> comJacobianPtr_;

  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  COMStates COMStates_;

  Eigen::VectorXd robotJointVelocity_;

  void updateCOMVelBounds_();
 
  /*
  double comVelXUpperBound_ = 0.0;
  double comVelXLowerBound_ = 0.0;

  double comVelYUpperBound_ = 0.0;
  double comVelYLowerBound_ = 0.0;
  */

  double pGain_ = 0.0;
  double dGain_ = 0.0;

  double poleOne_ = 0.0;
  double poleTwo_ = 0.0;

  double zUpperBound_ = 0.10;
  double zLowerBound_ = -0.07;
}; // End of struct name: ImpactAwareCOMVelConstraint

} // end of namespace mc_impact
