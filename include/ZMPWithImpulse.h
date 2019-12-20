#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include <McDynamicStability/McZMPArea.h>

#include "ConstraintUtils.h"

namespace mc_impact
{

template<typename Point>
struct ZMPWithImpulse : public mc_solver::InequalityConstraintRobot
{

  ///< ZMP defined with a set of points
  ZMPWithImpulse(mi_qpEstimator & predictor,
                 std::shared_ptr<mc_impact::McZMPArea<Point>> mcZMPAreaPtr,
                 const ImpactAwareConstraintParams<Point> & params);

  /*! \brief returns the upper bound of the amount of constraints for the QP to reserve memory accordingly.
   */
  inline int maxInEq() const override
  {
    if(updateMcZMPArea())
    {
      return getMcZMPArea()->getMaxNumVertex();
    }
    else
    {
      return static_cast<int>(getParams().zmpAreaVertexSet.size());
    }
  }

  inline std::string nameInEq() const override
  {
    return "ZMPWithImpulse";
  }

  inline const Eigen::MatrixXd & A() const override
  {
    return A_;
  }
  inline const Eigen::VectorXd & bInEq() const override
  {
    return b_;
  }
  /*!
   * \brief returns  the number of rows in realtime such that the QP can adjust the constraint size.
   */
  inline int nrInEq() const override
  {
    if(updateMcZMPArea())
    {
      return getMcZMPArea()->getNumVertex();
    }
    else
    {
      return static_cast<int>(getParams().zmpAreaVertexSet.size());
    }
  }
  void compute() override;

  void getInertialItems(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench);
  inline const Eigen::MatrixXd & getA()
  {
    return A_;
  }
  inline const Eigen::VectorXd & getb()
  {
    return b_;
  }
  inline const Eigen::MatrixXd & getZMP()
  {
    return A_zmp_;
  }

  inline bool updateMcZMPArea() const
  {
    return getParams().updateMcZMPArea;
  }
  inline const Eigen::Vector3d & getZMP_sensor()
  {
    if(getParams().debug)
      return zmpSensor_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_perturbation()
  {
    if(getParams().debug)
      return zmpPerturbation_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_prediction_feetforce()
  {
    if(getParams().debug)
      return zmpPrediction_feetforce_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_prediction_allforce()
  {
    if(getParams().debug)
      return zmpPrediction_allforce_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::VectorXd & getZMP_constraint_difference()
  {
    if(getParams().debug)
      return difference_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  /*
    inline const ZMPArea & getZMPArea(){
      return area_;
    }
    */
  inline const std::vector<Point> & getVertices()
  {
  if(updateMcZMPArea())
  {
    return getMcZMPArea()->getPolygonVertices();
  }else{
    return getParams().zmpAreaVertexSet;
  }
  }
  // Debugging:
  Point centeroid_;
  //Eigen::VectorXd slopeVec_ = Eigen::Vector3d::Zero();

  bool zmpTest_ = false;
  // Eigen::MatrixXd G_zmp_;
  // Eigen::VectorXd h_zmp_;

  
  inline void setIeqBlocks(const IeqConstraintBlocks & input)
  {
    ieqConstraintBlocks_ = input;
  }
  
  inline const IeqConstraintBlocks & getIeqBlocks() const
  {
    return ieqConstraintBlocks_;
  }

  Eigen::MatrixXd A_zmp_;
  bool pointInsideSupportPolygon(const Point & input);

  inline const std::shared_ptr<mc_impact::McZMPArea<Point>> getMcZMPArea() const
  {
    return mcZMPAreaPtr_;
  }

  inline const ImpactAwareConstraintParams<Point> & getParams() const
  {
    return params_;
  }

private:
  // Predictor
  mi_qpEstimator & predictor_;

  // Multi-contact ZMP area calculator:
  std::shared_ptr<mc_impact::McZMPArea<Point>> mcZMPAreaPtr_;

  // Timestep
  // double dt_;
  // Impact duration
  // double impact_dt_;

  // Alpha vector
  Eigen::VectorXd alpha_;

  // Name of the body of interest.
  // std::string bName_;
  // Name of the force sensor
  // std::string sName_;
  // dt * J_deltatau / impact_duration
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  IeqConstraintBlocks ieqConstraintBlocks_;

  //const std::vector<Point> iniVertexSet_;

  ImpactAwareConstraintParams<Point> params_;

  // ZMPArea area_;
  // bool allForce_;

  void calcZMP_();
  void computeMcZMPArea_();

  // double lowerSlope_;
  // double upperSlope_;
  // bool debug_;

  Eigen::Vector3d zmpSensor_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPerturbation_ = Eigen::Vector3d::Zero();

  Eigen::Vector3d zmpPrediction_allforce_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPrediction_feetforce_ = Eigen::Vector3d::Zero();

  Eigen::VectorXd difference_;
};

} // namespace mc_impact
