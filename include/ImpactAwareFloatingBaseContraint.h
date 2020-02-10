# pragma once 

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>


#include <McDynamicStability/McZMPArea.h>
#include <McDynamicStability/McComArea.h>
#include <McDynamicStability/McDCMArea.h>

#include "ConstraintUtils.h"


namespace mc_impact
{
struct ImpactAwareFloatingBaseConstraint: public mc_solver::InequalityConstraintRobot
{
  ///< ZMP defined with a set of points
  ImpactAwareFloatingBaseConstraint(mi_qpEstimator & predictor,
                 const mc_rbdyn::Robot & realRobot,
                 std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr,
                 std::shared_ptr<mc_impact::McComArea> mcComAreaPtr,
                 const ImpactAwareConstraintParams<Eigen::Vector2d> & params);




  /*! \brief We use this status to check what floating-base state do we constrain:
   * 1: only ZMP
   * 3: only DCM
   * 4: ZMP + DCM
   */
  int getConstrainingStatus() const{
    return constrainingStatus_; 
  }

  /*! \brief returns the upper bound of the amount of constraints for the QP to reserve memory accordingly.
  */
  inline int maxInEq() const override
  {
    switch(constrainingStatus_)
    {
     case 1: // Only constraining ZMP 
       return getMcZMPArea()->getMaxNumVertex();
     case 3: // Only constraining DCM
       return getMcDCMArea()->getMaxNumVertex();
     case 4: // Constraining both ZMP and DCM
       return (getMcZMPArea()->getMaxNumVertex() + getMcDCMArea()->getMaxNumVertex());
     default:
       throw std::runtime_error("The constrainingStatus is not set.");
    }
  }
  /*!
   * \brief returns  the number of rows in realtime such that the QP can adjust the constraint size.
   */
  inline int nrInEq() const override
  {
  switch(constrainingStatus_)
    {
     case 1: // Only constraining ZMP 
       return getMcZMPArea()->getNumVertex();
     case 3: // Only constraining DCM
       return getMcDCMArea()->getNumVertex();
     case 4: // Constraining both ZMP and DCM
       return (getMcZMPArea()->getNumVertex() + getMcDCMArea()->getNumVertex());
     default:
       throw std::runtime_error("The constrainingStatus is not set.");
    }
  }

  inline std::string nameInEq() const override
  {
    return "ImpactAwareFloatingBaseConstraint";
  }
  inline const Eigen::MatrixXd & A() const override
  {
    return A_;
  }
  inline const Eigen::VectorXd & bInEq() const override
  {
    return b_;
  }

  void compute() override;

  inline void setIeqBlocksZMP(const IeqConstraintBlocks & input)
  {
    if(getParams().constrainingZMP)
       ieqConstraintBlocksZMP_ = input;
    else
       throw std::runtime_error("ZMP inequality constraint is not set.");
  }

  inline const IeqConstraintBlocks & getIeqBlocksZMP() const
  {
     if(getParams().constrainingZMP)
       return ieqConstraintBlocksZMP_;
    else
       throw std::runtime_error("ZMP inequality constraint is not set.");

  }
   inline void setIeqBlocksDCM(const IeqConstraintBlocks & input)
  {
    if(getParams().constrainingDCM)
       ieqConstraintBlocksDCM_= input;
    else
       throw std::runtime_error("DCM inequality constraint is not set.");
  }

  inline const IeqConstraintBlocks & getIeqBlocksDCM() const
  {
     if(getParams().constrainingDCM)
       return ieqConstraintBlocksDCM_;
    else
       throw std::runtime_error("DCM inequality constraint is not set.");

  }

  inline const std::shared_ptr<const mc_impact::McZMPArea<Eigen::Vector2d>> getMcZMPArea() const
  {
    return mcZMPAreaPtr_;
  }

  inline const std::shared_ptr<const mc_impact::McComArea> getMcComArea() const
  {
    return mcComAreaPtr_;
  }

  inline const std::shared_ptr<const mc_impact::McDCMArea> getMcDCMArea() const
  {
    return mcDCMAreaPtr_;
  }

  inline const ImpactAwareConstraintParams<Eigen::Vector2d> & getParams() const
  {
    return params_;
  }
  inline const FloatingBaseStates & getFloatingBaseStates()const
  {
    return floatingBaseStates_; 
  } 
  inline const mc_rbdyn::Robot & realRobot() const
  {
   return realRobot_; 
  }
  /*! \brief get \omega = sqrt(9.8/com-z)
   */
  inline const double & getOmega()
  {
    return omega_;
  }
private:
  const mc_rbdyn::Robot & realRobot_;

  inline int  dof_() const
  {
   return realRobot_.mb().nrDof(); 
  }

  inline void calcOmega(const double & c_z)
  {
    assert(c_z > 0.0);
    omega_ = sqrt(9.8 / c_z);
  }

  double omega_;

  int getDCMRowNr_() const
  {
   switch(constrainingStatus_)
    {
     case 1: // Only constraining ZMP 
       throw std::runtime_error("Akding for DCM row number when only the ZMP constraint is used.");
     case 3: // Only constraining DCM
       return 0; 
     case 4: // Constraining both ZMP and DCM
       return getMcZMPArea()->getNumVertex();
     default:
       throw std::runtime_error("The constrainingStatus is not set.");
    }

  }

  // Predictor
  mi_qpEstimator & predictor_;
  // Multi-contact ZMP area calculator:
  std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr_;

  // Multi-contact Com area calculator:
  std::shared_ptr<mc_impact::McComArea> mcComAreaPtr_;

  // Multi-contact DCM area calculator:
  std::shared_ptr<mc_impact::McDCMArea> mcDCMAreaPtr_;

  void updateMcZMPAreas_(double height);
  void updateMcDCMAreas_();
  void getZMPBlocks(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench);

  int constrainingStatus_;

  void updateZMPConstraint_();
  void updateDCMConstraint_();

  Eigen::VectorXd robotJointVelocity_;

  // Left hand side of the ZMP and DCM constraint 
  Eigen::MatrixXd A_zmp_;
  Eigen::MatrixXd A_dcm_;
  void reset_();

  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  FloatingBaseStates floatingBaseStates_;
  void updateFloatingBaseState_();
  void calculateZMP_();

  IeqConstraintBlocks ieqConstraintBlocksZMP_;
  IeqConstraintBlocks ieqConstraintBlocksDCM_;

  std::shared_ptr<rbd::CoMJacobian> comJacobianPtr_;

  ImpactAwareConstraintParams<Eigen::Vector2d> params_;

}; // End of struct name: ImpactAwareFloatingBaseConstraint

} // end of namespace mc_impact
