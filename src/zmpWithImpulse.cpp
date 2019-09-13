#include "zmpWithImpulse.h"

namespace mc_impact
{

template<typename supportContact, typename Point>
zmpWithImpulse<supportContact, Point>::zmpWithImpulse(mi_qpEstimator & predictor,
                 const std::vector<supportContact> & supports,
                 double dt,
                 double impact_dt,
                 const std::vector<Point> & vertexSet,
                 bool allforce,
		 double lowerSlope, 
		 double upperSlope,
                 bool debug):InequalityConstraint(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt), impact_dt_(impact_dt), supports_(supports), iniVertexSet_(vertexSet), allForce_(allforce), debug_(debug){

  int numVertex = static_cast<int>(iniVertexSet_.size());
  A_zmp_ = Eigen::MatrixXd::Zero(numVertex, 6);

  //Eigen::MatrixXd G_zmp;
  //Eigen::VectorXd h_zmp;
   
  pointsToInequalityMatrix<Point>(iniVertexSet_, G_zmp_, h_zmp_, centeroid_, slopeVec_, lowerSlope, upperSlope);
  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  A_zmp_.block(0, 0, numVertex, 1) = G_zmp_.block(0, 1, numVertex, 1);
  /// A(:,1) = -G_x
  A_zmp_.block(0, 1, numVertex, 1) = - G_zmp_.block(0, 0, numVertex, 1);
  /// A(:,5) = h 
  A_zmp_.block(0, 5, numVertex, 1) = -h_zmp_; 


  int nDof = predictor_.getSimRobot().mb().nrDof();

  alpha_.resize(nDof);
  A_.resize(numVertex, nDof);
  b_.resize(numVertex);

  // Place holders
  /*
  area_.max_x = 0;
  area_.max_y = 0;
  area_.min_x = 0;
  area_.min_y = 0;
  */
}

/*
template <typename supportContact, typename Point>
zmpWithImpulse<supportContact, Point>::zmpWithImpulse(mi_qpEstimator & predictor,
                               const std::vector<supportContact> & supports,
                               double dt,
                               double impact_dt,
                               const ZMPArea & area,
                               bool allforce,
                               bool debug)
: InequalityConstraint(predictor.getSimRobot().robotIndex()), predictor_(predictor), dt_(dt), impact_dt_(impact_dt),
  supports_(supports), area_(area), allForce_(allforce), debug_(debug)
{

  // bName_ = bodyName;
  // sName_ = sensorName;
  A_zmp_ = Eigen::MatrixXd::Zero(4, 6);

  // - cy - max_x * fz
  A_zmp_(0, 1) = -1;
  A_zmp_(0, 5) = -area_.max_x;

  // cy + min_x * fz
  A_zmp_(1, 1) = 1;
  A_zmp_(1, 5) = area_.min_x;

  // cx - max_y * fz
  A_zmp_(2, 0) = 1;
  A_zmp_(2, 5) = -area_.max_y;

  // - cx + min_y * fz
  A_zmp_(3, 0) = -1;
  A_zmp_(3, 5) = area_.min_y;

  int nDof = predictor_.getSimRobot().mb().nrDof();

  alpha_.resize(nDof);
  // This needs double check
  A_.resize(4, nDof);
  b_.resize(4);

}
*/
template <typename supportContact, typename Point>
void zmpWithImpulse<supportContact, Point>::getInertialItems(Eigen::MatrixXd & sumJac, Eigen::Vector6d & exWrench)
{
  exWrench.setZero();
  int dof = predictor_.getSimRobot().mb().nrDof();
  sumJac.resize(6, dof);
  sumJac.setZero();
  // sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  // (1) Go through the bodies with contact

  for(auto idx = supports_.begin(); idx != supports_.end(); ++idx)
  {
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(idx->bodyName).inv();
    // exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().forceSensor(idx->sensorName).wrench().vector();
    exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(idx->bodyName).vector();

    //sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(idx->bodyName);

    
    sumJac.block(0, 0, 3, dof) += X_ee_0.dualMatrix().block(0, 3, 3, 3) * predictor_.getJacobianDeltaF(idx->bodyName);
    sumJac.block(3, 0, 3, dof) += X_ee_0.dualMatrix().block(3, 3, 3, 3) * predictor_.getJacobianDeltaF(idx->bodyName);
    //
  } // end of for

  // (2) Go through the impacts
  if(allForce_)
  {
    for(auto impactIdx = predictor_.getImpactModels().begin(); impactIdx != predictor_.getImpactModels().end();
        ++impactIdx)
    {
      sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();

      /*
      sumJac +=
          X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());

	  */
      
       sumJac.block(0, 0, 3, dof) +=
           X_ee_0.dualMatrix().block(0, 3, 3, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());

       sumJac.block(3, 0, 3, dof) +=
           X_ee_0.dualMatrix().block(3, 3, 3, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());
     

       // Add the hand force sensor measurement 
       //exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(impactIdx->second->getImpactBody()).vector();

    }
  }
}

template <typename supportContact, typename Point>
void zmpWithImpulse<supportContact, Point>::calcZMP_()
{
/*
  Eigen::VectorXd temp = (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alpha)
                          + rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD)
                                * predictor_.getImpactModels().begin()->second->getTimeStep());
				* */
  Eigen::VectorXd temp = predictor_.getImpactModels().begin()->second->getJointVel();
  double inv_t = (1 / predictor_.getImpactModels().begin()->second->getImpactDuration());

  Eigen::Vector6d local_exWrench;
  local_exWrench.setZero();
  int dof = predictor_.getSimRobot().mb().nrDof();
  Eigen::MatrixXd local_sumJac;
  local_sumJac.resize(6, dof);
  local_sumJac.setZero();
  // sva::PTransformd X_0_CoM = sva::PTransformd(predictor_.getRobot().com());
  // (1) Go through the bodies with contact

  for(auto idx = supports_.begin(); idx != supports_.end(); ++idx)
  {
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(idx->bodyName).inv();
    local_exWrench += X_ee_0.dualMatrix() * predictor_.getSimRobot().bodyWrench(idx->bodyName).vector();
    local_sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(idx->bodyName);

  } // end of for

  Eigen::Vector6d impulseForceSum;
  impulseForceSum.setZero();
  impulseForceSum = inv_t * local_sumJac * temp;
  double denominator = local_exWrench(5) + impulseForceSum(5);

  zmpPrediction_feetforce_.x() = -(local_exWrench(1) + impulseForceSum(1)) / denominator;
  zmpPrediction_feetforce_.y() = (local_exWrench(0) + impulseForceSum(0)) / denominator;

  // (2) Go through the impacts
  for(auto impactIdx = predictor_.getImpactModels().begin(); impactIdx != predictor_.getImpactModels().end();
      ++impactIdx)
  {
    sva::PTransformd X_ee_0 = predictor_.getSimRobot().bodyPosW(impactIdx->second->getImpactBody()).inv();

    local_sumJac += X_ee_0.dualMatrix().block(0, 3, 6, 3) * predictor_.getJacobianDeltaF(impactIdx->second->getImpactBody());
  }
  // Recalculate impulseSum:
  impulseForceSum = inv_t * local_sumJac * temp;
  denominator = local_exWrench(5) + impulseForceSum(5);

  zmpSensor_.x() = -(local_exWrench(1)) / denominator;
  zmpSensor_.y() = (local_exWrench(0)) / denominator;

  zmpPerturbation_.x() = -(impulseForceSum(1)) / denominator;
  zmpPerturbation_.y() = (impulseForceSum(0)) / denominator;

  zmpPrediction_allforce_ = zmpSensor_ + zmpPerturbation_;
}

template <typename supportContact, typename Point>
bool zmpWithImpulse<supportContact, Point>::pointInsideSupportPolygon(const Point & input){

  Point result = G_zmp_*input - h_zmp_;

  if (result[0] > 0)
	  return false;
  if (result[1] > 0)
	  return false;
  
  return true;
}

template <typename supportContact, typename Point>
void zmpWithImpulse<supportContact, Point>::computeAb()
{
  const auto & robot = predictor_.getSimRobot();
  Eigen::MatrixXd sumJac;
  Eigen::Vector6d sumWrench;
  getInertialItems(sumJac, sumWrench);


  A_ = (dt_ / impact_dt_) * A_zmp_ * sumJac;
  /*
    alpha_ =
          (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alpha));
    */
  rbd::paramToVector(robot.mbc().alpha, alpha_);
  b_ = -(A_zmp_ * sumWrench + A_zmp_ * sumJac * alpha_ / impact_dt_);

  if(debug_)
  {
    calcZMP_();
    difference_ = A_ * rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD) - b_;
  }
}


//The explicit instantiation
template struct mc_impact::zmpWithImpulse<zmpSupportContact, Eigen::Vector3d>;
template struct mc_impact::zmpWithImpulse<zmpSupportContact, Eigen::Vector2d>;


} // namespace mc_impact

