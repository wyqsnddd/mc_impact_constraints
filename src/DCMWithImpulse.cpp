#include "DCMWithImpulse.h"

namespace mc_impact
{

template<typename Point>
DCMWithImpulse<Point>::DCMWithImpulse(mi_qpEstimator & predictor,
                                      const mc_rbdyn::Robot & realRobot,
                                      std::shared_ptr<mc_impact::McZMPArea<Point>> mcZMPAreaPtr,
                                      std::shared_ptr<mc_impact::McComArea> mcComAreaPtr,
                                      const ImpactAwareConstraintParams<Point> & params)
: mc_solver::InequalityConstraintRobot(predictor.getSimRobot().robotIndex()), predictor_(predictor),
  mcZMPAreaPtr_(mcZMPAreaPtr), mcComAreaPtr_(mcComAreaPtr), realRobot_(realRobot), params_(params)
{

  int numVertex = static_cast<int>(getParams().zmpAreaVertexSet.size());
  /*
  pointsToInequalityMatrix<Point>(getParams().dcmAreaVertexSet, G_dcm_, h_dcm_, getParams().lowerSlope,
                                  getParams().upperSlope);
                                  */

  // A_dcm_.block(0, 0, numVertex, 1) = G_dcm_.block(0, 1, numVertex, 1);
  /// A(:,1) = -G_x
  // A_dcm_.block(0, 1, numVertex, 1) = -G_dcm_.block(0, 0, numVertex, 1);
  /// A(:,5) = h
  // A_dcm_.block(0, 5, numVertex, 1) = -h_dcm_;

  int nDof = predictor_.getSimRobot().mb().nrDof();
  alpha_.resize(nDof);
  A_.resize(numVertex, nDof);
  b_.resize(numVertex);

  difference_.resize(numVertex);
  difference_.setZero();

  comJacobianPtr_ = std::make_shared<rbd::CoMJacobian>(predictor_.getSimRobot().mb());

  calcOmega(predictor_.getSimRobot().com().z());

  // Create the Multi-contact DCM Area
  mcDCMAreaPtr_ = std::make_shared<mc_impact::McDCMArea>(mcZMPAreaPtr_, mcComAreaPtr_);

  // A_dcm_ = Eigen::MatrixXd::Zero(numVertex, 6);

  // Modify the ieqConstraintBlocks.
  pointsToInequalityMatrix(getParams().dcmAreaVertexSet, ieqConstraintBlocks_.G, ieqConstraintBlocks_.h,
                           getParams().lowerSlope, getParams().upperSlope);
  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  // A_dcm_.block(0, 0, numVertex, 1) = getIeqBlocks().G.block(0, 1, numVertex, 1);
  /// A(:,1) = -G_x
  // A_dcm_.block(0, 1, numVertex, 1) = -getIeqBlocks().G.block(0, 0, numVertex, 1);
  /// A(:,5) = h
  // A_dcm_.block(0, 5, numVertex, 1) = -getIeqBlocks().h;
}

template<typename Point>
bool DCMWithImpulse<Point>::pointInsideSupportPolygon(const Point & input)
{

  Eigen::VectorXd result = ieqConstraintBlocks_.G * input - ieqConstraintBlocks_.h;

  for(int ii = 0; ii < static_cast<int>(getParams().dcmAreaVertexSet.size()); ii++)
  {
    if(result(ii) > 0) return false;
  }

  return true;
}

template<typename Point>
void DCMWithImpulse<Point>::compute()
{

  mcZMPAreaPtr_->updateMcZMPArea(2.0);
  mcComAreaPtr_->updateMcComArea();
  // update the Multi-contact DCM area
  mcDCMAreaPtr_->updateMcDCMArea();

  // int numVertex = static_cast<int>(iniVertexSet_.size());
  // int numVertex = getMcDCMArea()->getNumVertex();

  // A_= Eigen::MatrixXd::Zero(numVertex, 6);

  // Set the inequality matrix blocks
  setIeqBlocks(getMcDCMArea()->getIeqConstraint());

  /// Needs to be checked carefully, compare to the 4 dim case
  // A(:,0) = G_y
  // A_dcm_.block(0, 0, numVertex, 1) = getIeqBlocks().G.block(0, 1, numVertex, 1);
  /// A(:,1) = -G_x
  // A_dcm_.block(0, 1, numVertex, 1) = -getIeqBlocks().G.block(0, 0, numVertex, 1);
  /// A(:,5) = h
  // A_dcm_.block(0, 5, numVertex, 1) = -getIeqBlocks().h;

  // const auto & robot = predictor_.getSimRobot();
  const auto & robot = realRobot_;
  int dof = predictor_.getSimRobot().mb().nrDof();
  Eigen::MatrixXd comJacobian = comJacobianPtr_->jacobian(robot.mb(), robot.mbc());
  // std::cout<<"comJacobian size is: "<<comJacobian.rows() << ", "<<comJacobian.cols()<<std::endl;
  Eigen::MatrixXd jacDcm = (comJacobian * predictor_.getJacobianDeltaAlpha()).block(0, 0, 2, dof);

  A_ = getParams().dt / getOmega() * getIeqBlocks().G * jacDcm;

  rbd::paramToVector(robot.mbc().alpha, alpha_);

  Eigen::Vector3d Com = robot.com();
  // Eigen::Vector3d
  ComVel_ = robot.comVelocity();

  dcm_ = Com.segment(0, 2) + ComVel_.segment(0, 2) / getOmega();

  predicted_dcm_jump_ = comJacobian.block(0, 0, 2, dof) * predictor_.getJointVelJump() / getOmega();

  predicted_dcm_ = dcm_ + predicted_dcm_jump_;

  b_ = (getIeqBlocks().h - getIeqBlocks().G * dcm_) - getIeqBlocks().G * jacDcm * alpha_ / getOmega();
  /*
   getOmega()*h_dcm_
   - G_dcm_*( dcm_*getOmega() + jacDcm*alpha_);

   */
  if(getParams().debug)
  {
    Eigen::VectorXd temp = (rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alpha)
                            + rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD)
                                  * predictor_.getImpactModels().begin()->second->getTimeStep());

    difference_ = A_ * rbd::dofToVector(predictor_.getSimRobot().mb(), predictor_.getSimRobot().mbc().alphaD) - b_;

    predictedComVelJump_ = comJacobian * predictor_.getJointVelJump();

    // The delta dq comparison:
    /*
     std::cout<<"The predictor delta dq: " <<std::endl<<predictor_.getJointVelJump().transpose()<<std::endl;

     std::cout<<"The calculated delta dq: "
     <<std::endl<<(predictor_.getJacobianDeltaAlpha()*temp).transpose()<<std::endl;

     std::cout<<"The difference is : " <<std::endl<<(predictor_.getJointVelJump() -
     predictor_.getJacobianDeltaAlpha()*temp).transpose()<<std::endl;
 */

    /*
        std::cout<<"The dcm constraint difference is: "<<difference_<<std::endl;
        std::cout<<"The dcm dq difference is: "<< G_dcm_*(dcm_ +
            comJacobian.block(0, 0, 2, dof)*predictor_.getJointVelJump() /getOmega()) - h_dcm_<<std::endl;

        std::cout<<"The dcm acc difference is: "<<G_dcm_*(dcm_  + jacDcm*temp/getOmega() ) - h_dcm_<<std::endl;


        std::cout<<"The dcm test is: "<<
          G_dcm_*dcm_  - h_dcm_<<std::endl;
        std::cout<<"The dcm jump test is: "<<
          G_dcm_*(dcm_ + predicted_dcm_jump_)  - h_dcm_<<std::endl;
        std::cout<<"The dcm prediction test is: "<<
          G_dcm_*(predicted_dcm_)  - h_dcm_<<std::endl;
    */
  }
}

// The explicit instantiation
// template struct mc_impact::DCMWithImpulse<Eigen::Vector3d>;
template struct mc_impact::DCMWithImpulse<Eigen::Vector2d>;

} // namespace mc_impact
