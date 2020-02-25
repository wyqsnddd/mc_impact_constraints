#pragma once

#include <mc_control/fsm/Controller.h>
#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ConstraintUtils.h"
#include <McDynamicStability/McComArea.h>
#include <McDynamicStability/McDCMArea.h>
#include <McDynamicStability/McZMPArea.h>

namespace mc_impact
{
struct ImpactAwareFloatingBaseConstraint : public mc_solver::InequalityConstraintRobot
{
  ///< ZMP defined with a set of points
  ImpactAwareFloatingBaseConstraint(
		  std::shared_ptr<mi_qpEstimator > predictorPtr,
		  std::shared_ptr<McContactSet> contactSetPtr,
		  const ImpactAwareConstraintParams<Eigen::Vector2d> & params);
  ~ImpactAwareFloatingBaseConstraint();


  /*! \brief We use this status to check what floating-base state do we constrain:
   * 1: only ZMP
   * 3: only DCM
   * 4: ZMP + DCM
   */
  int getConstrainingStatus() const
  {
    return constrainingStatus_;
  }

  inline bool zmpConstraintEnabled() const
  {

    if(getConstrainingStatus() == 3)
    {
      // Only DCM
      return false;
    }
    else
    {
      return true;
    }
  }

  inline bool dcmConstraintEnabled() const
  {
    if(getConstrainingStatus() == 1)
    {
      // Only ZMP
      return false;
    }
    else
    {
      // DCM, DCM + ZMP
      return true;
    }
  }

  /*! \brief returns the upper bound of the amount of constraints for the QP to reserve memory accordingly.
   */
  inline int maxInEq() const override
  {
    switch(constrainingStatus_)
    {
      case 1: // Only constraining ZMP
        return getZMPConstraintMaxSize_();
      case 3: // Only constraining DCM
        return getDCMConstraintMaxSize_();
      case 4: // Constraining both ZMP and DCM
        return (getZMPConstraintMaxSize_() + getDCMConstraintMaxSize_());
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
        return getZMPConstraintSize_();
      case 3: // Only constraining DCM
        return getDCMConstraintSize_();
      case 4: // Constraining both ZMP and DCM
        return (getZMPConstraintSize_() + getDCMConstraintSize_());
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
      ieqConstraintBlocksDCM_ = input;
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

  bool enabledMcDCMArea() const
  {
    return (getParams().constrainingDCM and getParams().updateMcDCMArea);
  }

  bool enabledMcZMPArea() const
  {
    return (getParams().updateMcZMPArea or getParams().updateMcDCMArea);
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

  inline const std::shared_ptr<const mc_impact::McComArea> getMcComArea() const
  {
    // McComArea will be updated if the McDCMArea is updated.
    if(getParams().updateMcDCMArea)
    {
      return mcComAreaPtr_;
    }
    else
    {
      LOG_ERROR_AND_THROW(std::runtime_error, "Asking for McDCMArea, which  is not initialized and updated.");
    }
  }

  inline const std::shared_ptr<const mc_impact::McDCMArea> getMcDCMArea() const
  {
    if(getParams().updateMcDCMArea)
    {
      return mcDCMAreaPtr_;
    }
    else
    {
      LOG_ERROR_AND_THROW(std::runtime_error, "Asking for McDCMArea, which  is not initialized and updated.");
    }
  }

  inline const ImpactAwareConstraintParams<Eigen::Vector2d> & getParams() const
  {
    return params_;
  }

  inline const FloatingBaseStates & getFloatingBaseStates() const
  {
    return floatingBaseStates_;
  }

  /*! \return reference to the robot (either pyhsices-engine simulated or real)
   */
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

  /*! \return if the 'samplePoint' is inside the multi-contact ZMP area.
   */
  bool pointInsideMcZMPArea(const Eigen::Vector3d & samplePoint) const;

  /*! \return if the 'samplePoint' is inside the multi-contact DCM area.
   */
  bool pointInsideMcDCMArea(const Eigen::Vector3d & samplePoint) const;

  void printDCMConstraintMatrix() const;
  void printZMPConstraintMatrix() const;

  /*! \brief Visualize Multi-contact areas in rviz.
   */
  void addMcAreasGuiItems(mc_control::fsm::Controller & ctl) const;

  /*! \brief Visualize floating-base states: Com, DCM and ZMP
   */
  void addFloatingBaseGuiItems(mc_control::fsm::Controller & ctl) const;

  /*! \brief Visualize Multi-contact surfaces, vertices in rviz.
   */
  void addMcContactGuiItems(mc_control::fsm::Controller & ctl) const;

  /*! \brif Add the logs.
   */
  void logFloatingBaseStates(mc_control::fsm::Controller & ctl);

  void removeLogFloatingBaseStates();

  inline std::shared_ptr<mi_qpEstimator> getPredictor()
  {
    return predictorPtr_; 
  }

private:
  // Predictor
  //mi_qpEstimator & predictor_;
  std::shared_ptr<mi_qpEstimator> predictorPtr_;
  ImpactAwareConstraintParams<Eigen::Vector2d> params_;
  std::shared_ptr<McContactSet> contactSetPtr_;

  inline void setLogger_()
  {
    loggersAdded_ = true;
  }
  inline bool checkLoggersAdded_() const
  {
    return loggersAdded_;
  }
  bool loggersAdded_ = false;

  inline mc_control::fsm::Controller * getControllerWithLogger_() 
  {
    if(controllerWithLoggerPtr_ != nullptr)
    {
      return controllerWithLoggerPtr_;
    }
    else
    {
      throw std::runtime_error("Getting controllerWithLoggerPtr_ without setting."); 
    }
  }

  inline void setControllerWithLogger_(mc_control::fsm::Controller & ctl) 
  {
    controllerWithLoggerPtr_ = &ctl; 
  }

  mc_control::fsm::Controller *  controllerWithLoggerPtr_;


  const mc_rbdyn::Robot & robot_;

  inline int dof_() const
  {
    return robot().mb().nrDof();
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
        throw std::runtime_error("Asking for DCM row number when only the ZMP constraint is used.");
      case 3: // Only constraining DCM
        return 0;
      case 4: // Constraining both ZMP and DCM
        if(getParams().updateMcZMPArea)
          return getMcZMPArea()->getNumVertex();
        else
          return static_cast<int>(getParams().zmpAreaVertexSet.size());
      default:
        throw std::runtime_error("The constrainingStatus is not set.");
    }
  }

  // Multi-contact ZMP area calculator:
  std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr_;

  // Multi-contact Com area calculator:
  std::shared_ptr<mc_impact::McComArea> mcComAreaPtr_;

  // Multi-contact DCM area calculator:
  std::shared_ptr<mc_impact::McDCMArea> mcDCMAreaPtr_;

  void updateMcAreas_(double height);
  void updateMcZMPAreas_(double height);
  void updateMcDCMAreas_();

  inline int getZMPConstraintSize_() const
  {
    if(getParams().updateMcZMPArea)
    {
      return getMcZMPArea()->getNumVertex();
    }
    else
    {
      return static_cast<int>(getParams().zmpAreaVertexSet.size());
    }
  }

  inline int getDCMConstraintSize_() const
  {
    if(getParams().updateMcDCMArea)
    {
      return getMcDCMArea()->getNumVertex();
    }
    else
    {
      return static_cast<int>(getParams().dcmAreaVertexSet.size());
    }
  }

  inline int getZMPConstraintMaxSize_() const
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

  inline int getDCMConstraintMaxSize_() const
  {
    if(getParams().updateMcDCMArea)
    {
      return getMcDCMArea()->getMaxNumVertex();
    }
    else
    {
      return static_cast<int>(getParams().dcmAreaVertexSet.size());
    }
  }
  void getZMPBlocks(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench);

  int constrainingStatus_;

  void updateZMPConstraint_();
  void updateDCMConstraint_();

  Eigen::VectorXd robotJointVelocity_;

  // Left hand side of the ZMP and DCM constraint
  Eigen::MatrixXd A_zmp_;
  Eigen::MatrixXd A_dcm_;
  void fixedSupportPolygonSetup_();

  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  // Needs to be reset in each iteration if we update the Multi-contact areas in realtime.
  void updateAbMatrixsize_();

  FloatingBaseStates floatingBaseStates_;
  void updateFloatingBaseState_();
  void calculateZMP_();

  IeqConstraintBlocks ieqConstraintBlocksZMP_;
  IeqConstraintBlocks ieqConstraintBlocksDCM_;

  std::shared_ptr<rbd::CoMJacobian> comJacobianPtr_;


}; // End of struct name: ImpactAwareFloatingBaseConstraint

} // end of namespace mc_impact
