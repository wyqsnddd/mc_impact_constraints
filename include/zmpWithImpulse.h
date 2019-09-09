#pragma once

#include <mc_prediction/mi_qpEstimator.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>
#include <mc_solver/QPSolver.h>

# include "constraintUtils.h"

namespace mc_impact
{
template<typename supportContact, typename Point>
struct zmpWithImpulse : public mc_solver::InequalityConstraint
{

  /// ZMP defined by a rectangle
  zmpWithImpulse(
		 mi_qpEstimator & predictor,
                 const std::vector<supportContact> & supports,
                 double dt,
                 double impact_dt,
                 const ZMPArea & area,
                 bool allforce = true,
                 bool debug = false);

 
  /// ZMP defined with a set of points
  zmpWithImpulse( 
		 mi_qpEstimator & predictor,
                 const std::vector<supportContact> & supports,
                 double dt,
                 double impact_dt,
                 const std::vector<Point> & vertexSet,
                 bool allforce = true,
                 bool debug = false);

  inline int maxInEq() const override
  {
    return static_cast<int>(iniVertexSet_.size());
  }

  inline std::string nameInEq() const override
  {
    return "zmpWithImpulse";
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
  inline const Eigen::Vector4d & getZMP_constraint_difference()
  {
    if(debug_)
      return difference_;
    else
      throw std::runtime_error("zmp constraint not in debug mode.");
  }

  inline const ZMPArea & getZMPArea(){
    return area_; 
  }
  inline const std::vector<Point> & getVertices(){
    return iniVertexSet_; 
  } 
// Debugging: 
  Point centeroid_;
  Eigen::VectorXd slopeVec_ = Eigen::Vector3d::Zero();

  bool zmpTest_ = false;
  Eigen::MatrixXd G_zmp_;
  Eigen::VectorXd h_zmp_;

  Eigen::MatrixXd A_zmp_;
  bool pointInsideSupportPolygon(const Point & input);
private:
  // Predictor
  mi_qpEstimator & predictor_;
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

  const std::vector<Point> iniVertexSet_;

  ZMPArea area_; 
  bool allForce_;
  bool debug_;

  void calcZMP_();
  Eigen::Vector3d zmpSensor_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPerturbation_ = Eigen::Vector3d::Zero();

  Eigen::Vector3d zmpPrediction_allforce_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d zmpPrediction_feetforce_ = Eigen::Vector3d::Zero();

  Eigen::Vector4d difference_ = Eigen::Vector4d::Zero();
};

} // namespace mc_impact
