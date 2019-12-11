#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "constraintUtils.h"

#include <mc_dynamicStability/mc_zmp_area.h>

namespace mc_impact
{




template<typename supportContact, typename Point>
struct ZMPWithImpulse : public mc_solver::InequalityConstraintRobot
{

  /// ZMP defined with a set of points
  ZMPWithImpulse(mi_qpEstimator & predictor,
		  std::shared_ptr<mc_impact::McZMPArea<Point> > mcZMPAreaPtr,
                 const std::vector<supportContact> & supports,
                 double dt,
                 double impact_dt,
                 const std::vector<Point> & vertexSet,
                 bool allforce = true,
		 bool updateMcZMPArea = true,
                 double lowerSlope = 0.01,
                 double upperSlope = 100.0,
                 bool debug = false);

  /*! \brief returns the upper bound of the amount of constraints for the QP to reserve memory accordingly.
   */
  inline int maxInEq() const override
  {
    if(updateMcZMPArea())
    {
      return getMcZMPArea()->getMaxNumVertex(); 
    }else{
      return static_cast<int>(iniVertexSet_.size());
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
    }else{
      return static_cast<int>(iniVertexSet_.size());
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
    return updateMcZMPArea_; 
  }
  inline const Eigen::Vector3d & getZMP_sensor()
  {
    if(debug_)
      return zmpSensor_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_perturbation()
  {
    if(debug_)
      return zmpPerturbation_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_prediction_feetforce()
  {
    if(debug_)
      return zmpPrediction_feetforce_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::Vector3d & getZMP_prediction_allforce()
  {
    if(debug_)
      return zmpPrediction_allforce_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }
  inline const Eigen::VectorXd & getZMP_constraint_difference()
  {
    if(debug_)
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
    return iniVertexSet_;
  }
  // Debugging:
  Point centeroid_;
  Eigen::VectorXd slopeVec_ = Eigen::Vector3d::Zero();

  bool zmpTest_ = false;
  //Eigen::MatrixXd G_zmp_;
  //Eigen::VectorXd h_zmp_;

  inline void setIeqBlocks(const ieqConstraintBlocks & input)
  {
    ieqConstraintBlocks_ = input; 
  }
  inline const ieqConstraintBlocks & getIeqBlocks() const
  {
    return ieqConstraintBlocks_; 
  }

  Eigen::MatrixXd A_zmp_;
  bool pointInsideSupportPolygon(const Point & input);

  inline const std::shared_ptr<mc_impact::McZMPArea<Point> > getMcZMPArea() const
  {
    return mcZMPAreaPtr_; 
  }
private:
  // Predictor
  mi_qpEstimator & predictor_;

  // Multi-contact ZMP area calculator: 
  std::shared_ptr<mc_impact::McZMPArea<Point> > mcZMPAreaPtr_;

  // Timestep
  double dt_;
  // Impact duration
  double impact_dt_;
  // Alpha vector
  Eigen::VectorXd alpha_;
  std::vector<supportContact> supports_;
  // Name of the body of interest.
  // std::string bName_;
  // Name of the force sensor
  // std::string sName_;
  // dt * J_deltatau / impact_duration
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  ieqConstraintBlocks  ieqConstraintBlocks_;

  const std::vector<Point> iniVertexSet_;

  // ZMPArea area_;
  bool allForce_;

  bool updateMcZMPArea_;
  void calcZMP_();
  void computeMcZMPArea_(); 


  double lowerSlope_;
  double upperSlope_;
  bool debug_;


    Eigen::Vector3d zmpSensor_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPerturbation_ = Eigen::Vector3d::Zero();

  Eigen::Vector3d zmpPrediction_allforce_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPrediction_feetforce_ = Eigen::Vector3d::Zero();

  Eigen::VectorXd difference_; 

};

} // namespace mc_impact
