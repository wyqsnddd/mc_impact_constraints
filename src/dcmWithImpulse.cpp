#include "dcmWithImpulse.h"

namespace mc_impact {


template<typename supportContact, typename Point>
dcmWithImpulse<supportContact, Point>::dcmWithImpulse(
		const mc_rbdyn::Robot & realRobot,
		mi_qpEstimator & predictor, 
                 const std::vector<supportContact> & supports,
                 double dt,
                 double impact_dt,
                 const std::vector<Point> & vertexSet,
		 double lowerSlope, 
		 double upperSlope,
		 bool debug
                 ):InequalityConstraint(predictor.getSimRobot().robotIndex()), realRobot_(realRobot), predictor_(predictor), dt_(dt), impact_dt_(impact_dt), supports_(supports), iniVertexSet_(vertexSet), debug_(debug) {

  int numVertex = static_cast<int>(iniVertexSet_.size());
  A_dcm_ = Eigen::MatrixXd::Zero(numVertex, 6);

  pointsToInequalityMatrix<Point>(iniVertexSet_, G_dcm_, h_dcm_, centeroid_, slopeVec_, lowerSlope, upperSlope);
 
  A_dcm_.block(0, 0, numVertex, 1) = G_dcm_.block(0, 1, numVertex, 1);
  /// A(:,1) = -G_x
  A_dcm_.block(0, 1, numVertex, 1) = - G_dcm_.block(0, 0, numVertex, 1);
  /// A(:,5) = h 
  A_dcm_.block(0, 5, numVertex, 1) = -h_dcm_; 


  int nDof = predictor_.getSimRobot().mb().nrDof();
  alpha_.resize(nDof);
  A_.resize(numVertex, nDof);
  b_.resize(numVertex);

  difference_.resize(numVertex);
  difference_.setZero();

  comJacobianPtr_ = std::make_shared<rbd::CoMJacobian>(predictor_.getSimRobot().mb());

  calcOmega( predictor_.getSimRobot().com().z() );
}




template <typename supportContact, typename Point>
bool dcmWithImpulse<supportContact, Point>::pointInsideSupportPolygon(const Point & input){

  Eigen::VectorXd result = G_dcm_*input - h_dcm_;

  for (int ii = 0; ii <= static_cast<int>(iniVertexSet_.size()); ii++){
    if (result(ii) >0)
	    return false;
  }

  return true;
}

template <typename supportContact, typename Point>
void dcmWithImpulse<supportContact, Point>::computeAb()
{
  //const auto & robot = predictor_.getSimRobot();
  const auto & robot = realRobot_; 
  int dof = predictor_.getSimRobot().mb().nrDof();
  Eigen::MatrixXd comJacobian = comJacobianPtr_->jacobian(robot.mb(), robot.mbc());
  //std::cout<<"comJacobian size is: "<<comJacobian.rows() << ", "<<comJacobian.cols()<<std::endl;
  Eigen::MatrixXd jacDcm =   (comJacobian * predictor_.getJacobianDeltaAlpha()).block(0, 0, 2, dof);

  A_ = dt_/getOmega()  * G_dcm_ * jacDcm; 
  

  rbd::paramToVector(robot.mbc().alpha, alpha_);

  Eigen::Vector3d Com = robot.com();
  Eigen::Vector3d ComVel = robot.comVelocity();

  dcm_ =  Com.segment(0, 2) + ComVel.segment(0,2)/getOmega();

  predicted_dcm_jump_ =  comJacobian.block(0, 0, 2, dof)*predictor_.getJointVelJump() /getOmega(); 

  predicted_dcm_ = dcm_ + predicted_dcm_jump_; 

  b_ = (h_dcm_ - G_dcm_*dcm_)  - G_dcm_*jacDcm*alpha_/getOmega();
	 /* 
	  getOmega()*h_dcm_ 
	  - G_dcm_*( dcm_*getOmega() + jacDcm*alpha_);

	  */
  if(debug_)
  {
    Eigen::VectorXd temp = (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alpha)
                          + rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD)
                                * predictor_.getImpactModels().begin()->second->getTimeStep());

    difference_ = A_ * rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD) - b_;

    predictedComVelJump_ = comJacobian*predictor_.getJointVelJump(); 

   // The delta dq comparison: 
   /*
    std::cout<<"The predictor delta dq: " <<std::endl<<predictor_.getJointVelJump().transpose()<<std::endl;

    std::cout<<"The calculated delta dq: " <<std::endl<<(predictor_.getJacobianDeltaAlpha()*temp).transpose()<<std::endl;

    std::cout<<"The difference is : " <<std::endl<<(predictor_.getJointVelJump() - predictor_.getJacobianDeltaAlpha()*temp).transpose()<<std::endl;
*/


    std::cout<<"The dcm constraint difference is: "<<difference_<<std::endl;
    //std::cout<<"A_ size: "<<A_.rows()<<", "<<A_.cols()<<std::endl;
    //std::cout<<"b_ size: "<<b_.rows()<<", "<<b_.cols()<<std::endl;
    std::cout<<"The dcm dq difference is: "<< G_dcm_*(dcm_ + 
		    comJacobian.block(0, 0, 2, dof)*predictor_.getJointVelJump() /getOmega()) - h_dcm_<<std::endl;

    std::cout<<"The dcm acc difference is: "<<G_dcm_*(dcm_  + jacDcm*temp/getOmega() ) - h_dcm_<<std::endl;
    

    std::cout<<"The dcm test is: "<<
	    G_dcm_*dcm_  - h_dcm_<<std::endl;
    std::cout<<"The dcm jump test is: "<<
	    G_dcm_*(dcm_ + predicted_dcm_jump_)  - h_dcm_<<std::endl;
    std::cout<<"The dcm prediction test is: "<<
	    G_dcm_*(predicted_dcm_)  - h_dcm_<<std::endl;



  }

}

//The explicit instantiation
template struct mc_impact::dcmWithImpulse<zmpSupportContact, Eigen::Vector3d>;
template struct mc_impact::dcmWithImpulse<zmpSupportContact, Eigen::Vector2d>;



} // end of namespace 
