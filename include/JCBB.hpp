/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung, Felipe Inostroza
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 *       Universidad de Chile, nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef JCBB_HPP
#define JCBB_HPP

#include <boost/foreach.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/shared_ptr.hpp>
#include <cstddef>
#include <float.h>
#include <queue>
#include <vector>

#include "MeasurementModel.hpp"
#include "Tree.hpp"

namespace rfs{

  /**
   * \class JCBB_TreeNode
   * \brief Node for the JCBB interpretation tree
   * Node for the JCBB interpretation tree
   */
  class JCBB_TreeNode: public TreeNode<JCBB_TreeNode>
  {
  public:

    /**
     * Default constructor
     */
    JCBB_TreeNode();

    /**
     * Destructor
     */
    ~JCBB_TreeNode();

    /** 
     * Less than operator for use with priorty queue
     */
    bool operator < (JCBB_TreeNode &other) const;

    /**
     * Set the data association represented by this node
     * \param[in] z_idx measurement index
     * \param[in] m_idx landmark index
     * \param[in] dist_mahalanobis2 squared Mahalanobis distance for the associations
     * represented by this and all parent nodes
     */
    void setAssociation(int z_idx, int m_idx, double dist_mahalanobis2);

    /**
     * Get the association represented by this node
     * \param[out] z_idx measurement index
     * \param[out] m_idx landmark index
     */ 
    void getAssociation(int &z_idx, int &m_idx);

    /**
     * Get the joint Mahalanobis distance squared of all associations up to the root
     * \return Mahalanobis distance squared
     */
    double getMahalanobisDist2();

    /**
     * Get the measurement innovation of all associations up to the root
     * \return pointer to innovation vector
     */
    ::Eigen::VectorXd *getInnovation();

    /**
     * Get the measurement innovation covariance of all assoications up to the root
     * \return pointer to innovation covariance
     */
    ::Eigen::MatrixXd* getInnovationCovInv();

  private:

    int z_idx_; /**< measurement index */
    int m_idx_; /**< landmark index */
    double dist_mahalanobis2_; /**< joint Mahalanobis distance for all associations up to the root */
    ::Eigen::MatrixXd C_inverse_; /**< inverse innovation covariance for all associations up to the root */
    ::Eigen::VectorXd innov_; /**< measurement innovation for all associations up to the root */

  };

  /**
   * \class JCBB
   * \brief Joint compatibility branch and bound data association algorithm
   *
   * An implementation of the Joint compatibility branch and bound data association algorithm
   * Reference:  Neira, J. Tardos; Tardos, J.D. (2001), "Data association in stochastic mapping 
   * using the joint compatibility test", Robotics and Automation, IEEE Transactions on 17 (6): 
   * pp 890–897, doi:10.1109/70.976019.
   */
  template<class MeasurementModelType>
  class JCBB
  {

  public:

    typedef typename MeasurementModelType::TPose TPose;
    typedef typename MeasurementModelType::TLandmark TLandmark;
    typedef typename MeasurementModelType::TMeasurement TMeasurement;
    typedef typename MeasurementModelType::TJacobianPose TJacobianPose;
    typedef typename MeasurementModelType::TJacobianLmk TJacobianLmk;
    typedef boost::shared_ptr<JCBB_TreeNode> JCBB_TreeNodePtr;

    /**
     * Constructor: Setup and initiate the JCBB algorithm. For correlated robot and landmark estimates,
     * use the estCovDense parameter to set the estimate erro covariance matrix. If this 
     * is not set, it is assumed that the covariance matrix is block diagonal, and that
     * the individual robot and landmark covariances can be found in the TPose and 
     * TLandmark objects.
     * \param[in] confidenceItvl confidence interval for gating
     * \param[in] measurementModel pointer to the measurement model. The Jacobian of the
     * meausrement function w.r.t. the robot pose and landmark positions must be defined.
     * With the exception that if the robot pose estimate covariance is zero, than the 
     * Jacobian w.r.t the robot pose does not need to be defined.
     * \param[in] measurements pointer to a vector of measurements
     * \param[in] robot pointer to the robot / sensor pose.
     * \param[in] landmarks pointer to a vector of landmark pointers
     * \param[in] estCovDense the dense estimate error covariance matrix.
     */
    JCBB(double confidenceItvl,
	 MeasurementModelType *measurementModel,
	 std::vector< TMeasurement > *measurements,
	 TPose const * robot,
	 std::vector< TLandmark* > *landmarks,
	 ::Eigen::MatrixXd const *estCovDense = NULL);

    /**
     * Destructor
     */
    ~JCBB();

    /**
     * Get the landmark index associated with a measurement as determined by JCBB
     * \param[in] z_idx measurement index
     * \return the landmark index, -1 if unassociated
     */
    int getAssociation(int z_idx);


  private:

    double confidence_interval_; /**< Chi-square confidence interval for gating  */

    JCBB_TreeNodePtr iTreeRoot_; /**< Interpretation tree */
    
    MeasurementModelType *measurementModel_; /**< measurement model */
    std::vector<TMeasurement> *Z_actual_ptr_; /**< vector of actual measurements */
    std::vector<TMeasurement> Z_predict_; /**< vector expected measurements */
    std::vector< boost::shared_ptr< TJacobianPose > > jacobians_wrt_pose_; /**< Jacobian of of measurement model w.r.t. the robot pose */
    std::vector< boost::shared_ptr< TJacobianLmk  > > jacobians_wrt_lmk_; /**< Jacobian of the measurement model w.r.t. landmark positions */

    TPose const *robot_; /**< robot pose */ 
    std::vector<TLandmark*> *landmarks_; /**< vector of landmark positions (and covariance) */
    bool isEstCovBlockDiag_; /**< flag for block diagonal or dense estimate covariance matrix */
    const ::Eigen::MatrixXd *estCovDense_; /**< Dense estimate covariance matrix */

    std::priority_queue<JCBB_TreeNodePtr> branchList_; /**< priority queue of nodes in the interpretation tree to branch */

    int nAssociations_best_; /**< the number of assoications in the best hypothesis */
    double dist_mahalanobis2_best_; /**< the mahalanobis distance of the best hypothesis */
    JCBB_TreeNodePtr node_best_; /**< the node in the interpretation tree associated with the best hypothesis */
    std::vector<int> hypothesis_best_; /**< data association indices for the best hypothesis */

    /**
     * calculate expected measurements and covariance 
     */
    void calcExpectedMeasurements();
    
    /**
     * Branch a node of the interpretation tree
     * \param[in] node Node to be branched
     */
    void branch(JCBB_TreeNodePtr node);

  };

} // namespace rfs


// implementation JCBB_TreeNode

namespace rfs{

  JCBB_TreeNode::JCBB_TreeNode(){

  }

  JCBB_TreeNode::~JCBB_TreeNode(){

  }

  bool JCBB_TreeNode::operator<(JCBB_TreeNode &other) const{
    if(this->dist_mahalanobis2_ < other.dist_mahalanobis2_)
      return true;
    else
      return false;
  }

  void JCBB_TreeNode::setAssociation(int z_idx, int m_idx, double dist_mahalanobis2){
    z_idx_ = z_idx;
    m_idx_ = m_idx;
    dist_mahalanobis2_ = dist_mahalanobis2;
  }
  
  void JCBB_TreeNode::getAssociation(int &z_idx, int &m_idx){
    z_idx = z_idx_;
    m_idx = m_idx_;
  }

  double JCBB_TreeNode::getMahalanobisDist2(){
    return dist_mahalanobis2_;
  }

  ::Eigen::VectorXd* JCBB_TreeNode::getInnovation(){
    return &innov_;
  }

  ::Eigen::MatrixXd* JCBB_TreeNode::getInnovationCovInv(){
    return &C_inverse_;
  }

}


// implementation JCBB

namespace rfs{

  template<class MeasurementModelType>
  JCBB<MeasurementModelType>::JCBB(double confidenceItvl,
				   MeasurementModelType *measurementModel,
				   std::vector< TMeasurement > *measurements,
				   TPose const *robot,
				   std::vector< TLandmark* > *landmarks,
				   ::Eigen::MatrixXd const *estCovDense){

    
    iTreeRoot_ = JCBB_TreeNodePtr(new JCBB_TreeNode() );
    
    confidence_interval_ = confidenceItvl;
    measurementModel_ = measurementModel;
    Z_actual_ptr_ = measurements;
    robot_ = robot;
    landmarks_ = landmarks;
    isEstCovBlockDiag_ = true;
    if(estCovDense != NULL){
      isEstCovBlockDiag_ = false;
      estCovDense_ = estCovDense;
    }
    nAssociations_best_ = 0;
    dist_mahalanobis2_best_ = DBL_MAX;

    if(confidence_interval_ > 1)
      return;

    calcExpectedMeasurements();

    iTreeRoot_->setAssociation(-1, -1, 0);
    branchList_.push(iTreeRoot_);
    while(!branchList_.empty()){
      JCBB_TreeNodePtr topNode = branchList_.top();
      branchList_.pop();
      branch( topNode );
    }

    hypothesis_best_.clear();
    hypothesis_best_.resize( Z_actual_ptr_->size(), -1 );
    JCBB_TreeNodePtr n = node_best_;
    while(n != iTreeRoot_){
      int z_idx, m_idx;
      n->getAssociation(z_idx, m_idx);
      hypothesis_best_[z_idx] = m_idx;
      n = n->getParent();
    }

  }

  template<class MeasurementModelType>
  JCBB<MeasurementModelType>::~JCBB(){

  }

  template<class MeasurementModelType>
  int JCBB<MeasurementModelType>::getAssociation(int z_idx){
    return hypothesis_best_[z_idx];
  }

  template<class MeasurementModelType>
  void JCBB<MeasurementModelType>::calcExpectedMeasurements(){

    int nLandmarks = landmarks_->size();
    TMeasurement z_predict; 
    Z_predict_.clear();
    Z_predict_.reserve( nLandmarks );
    jacobians_wrt_pose_.reserve( nLandmarks );
    jacobians_wrt_lmk_.reserve( nLandmarks );

    BOOST_FOREACH( TLandmark* lmkPtr, *landmarks_){
      
      boost::shared_ptr< TJacobianPose > H_pose_ptr(new TJacobianPose); 
      boost::shared_ptr< TJacobianLmk  > H_lmk_ptr(new TJacobianLmk); 
      
      measurementModel_->measure(*robot_, *lmkPtr, z_predict, H_lmk_ptr.get(), H_pose_ptr.get());
      Z_predict_.push_back(z_predict);
      jacobians_wrt_pose_.push_back( H_pose_ptr );
      jacobians_wrt_lmk_.push_back( H_lmk_ptr );

    }
  }

  template<class MeasurementModelType>
  void JCBB<MeasurementModelType>::branch(JCBB_TreeNodePtr node){

    // Check if all measurements have been used up
    int z_idx_current;
    int m_idx;
    node->getAssociation(z_idx_current, m_idx);
    if(z_idx_current != -1 && z_idx_current + 1 >= Z_actual_ptr_->size()){
      return;
    }

    // Check previous associations
    int z_idx;
    JCBB_TreeNodePtr node_current = node;
    std::vector<int> isUsed_lmk_idx(landmarks_->size(), 0);
    std::vector<int> used_lmk_indices;
    while(node_current->getParent() != iTreeRoot_ && node_current->getParent() != NULL){
      node_current->getAssociation(z_idx, m_idx);
      if(m_idx >= 0){
	isUsed_lmk_idx[ m_idx ] = 1;
	used_lmk_indices.push_back( m_idx );
      }
      node_current = node_current->getParent();
    }
    int nAssociations = used_lmk_indices.size();
    ::Eigen::VectorXd *innov_m = node->getInnovation();
    double d2_m = node->getMahalanobisDist2();
    
    // Determine next measurement to associate
    node->getAssociation(z_idx, m_idx);
    z_idx++;

    // Branch (for unused landmark indices), if not rejected by gating
    m_idx = 0;;
    BOOST_FOREACH(int isUsed, isUsed_lmk_idx){
      if( isUsed == 0 ){ // m_idx has not been associated with another measurement
	
	// innovation
	typename TMeasurement::Vec z_actual = Z_actual_ptr_->at(z_idx).get();
	typename TMeasurement::Vec z_predict = Z_predict_[m_idx].get();
	typename TMeasurement::Vec innov = z_actual - z_predict; 
	typename TMeasurement::Mat Q;
	measurementModel_->getNoise(Q);

	int robotDim = robot_->getNDim();
	int lmkDim = landmarks_->at(0)->getNDim();
	int measureDim = z_actual.size();
	int nPreviousMeasurements = used_lmk_indices.size();
	double dist_mahalanobis2 = -1;

	::Eigen::MatrixXd W;
	W.resize(measureDim, measureDim * nPreviousMeasurements);

	::Eigen::Matrix<double, TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> C_current;

	bool isRobotPoseHasUncertainty = true;

	if(!isEstCovBlockDiag_){ // Estimate covariance is not block diagonal
	
	  ::Eigen::MatrixXd P_r_rows = estCovDense_->topRows(robotDim); // first row-block (robotDim rows) of estimate cov	
	  ::Eigen::MatrixXd P_lmk_rows = estCovDense_->block(robotDim + m_idx*lmkDim, 0, lmkDim, estCovDense_->cols() ); // the row-block of landmark m_idx 
	  ::Eigen::MatrixXd HP = *(jacobians_wrt_pose_[m_idx]) * P_r_rows + *(jacobians_wrt_lmk_[m_idx]) * P_lmk_rows;
	  ::Eigen::MatrixXd HP_cols_r = HP.leftCols(robotDim); // first col-block (robotDim cols) of HP

	  int W_block_idx = 0;
	  BOOST_REVERSE_FOREACH(int m_idx_used, used_lmk_indices){	 
	    ::Eigen::MatrixXd HP_cols_lmk = HP.block(0, robotDim + m_idx_used*lmkDim, HP.rows(), lmkDim);  // Pick off the col-block of HP for landmark m_idx_used
	    W.block(0, measureDim * W_block_idx, measureDim, measureDim) = HP_cols_r * jacobians_wrt_pose_[m_idx_used]->transpose() + 
	      HP_cols_lmk * jacobians_wrt_lmk_[m_idx_used]->transpose() + Q;
	    W_block_idx++;
	  }

	  C_current = *(jacobians_wrt_pose_[m_idx]) * estCovDense_->block(0, 0, robotDim, robotDim) * (jacobians_wrt_pose_[m_idx])->transpose() + 
	    *(jacobians_wrt_lmk_[m_idx]) * estCovDense_->block(robotDim + m_idx*lmkDim, robotDim + m_idx*lmkDim, lmkDim, lmkDim) * (jacobians_wrt_lmk_[m_idx])->transpose() + Q;

	}else{ // Estimate covariance is block diagonal

	  ::Eigen::Matrix<double, TPose::Vec::RowsAtCompileTime, TPose::Vec::RowsAtCompileTime> P_rr = robot_->getCov();
	  ::Eigen::Matrix<double, TLandmark::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime> P_mm = landmarks_->at(m_idx)->getCov();

	  if(!P_rr.isZero()){
	    for(int W_block_idx = 0; W_block_idx < nPreviousMeasurements; W_block_idx++){		
	      W.block(0, measureDim * W_block_idx, measureDim, measureDim) = *(jacobians_wrt_pose_[m_idx]) * P_rr * jacobians_wrt_pose_[m_idx]->transpose() + Q;
	    }

	    C_current = *(jacobians_wrt_pose_[m_idx]) * P_rr * (jacobians_wrt_pose_[m_idx])->transpose() + 
	                *(jacobians_wrt_lmk_[m_idx])  * P_mm * (jacobians_wrt_lmk_[m_idx])->transpose() + Q;

	  }else{
	    W.setZero();

	    C_current = *(jacobians_wrt_lmk_[m_idx]) * P_mm * (jacobians_wrt_lmk_[m_idx])->transpose() + Q;

	    isRobotPoseHasUncertainty = false;
	  }

	}

	::Eigen::Matrix<double, TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> N;
	::Eigen::MatrixXd L;
	::Eigen::MatrixXd K;

	if(nPreviousMeasurements != 0 && (isRobotPoseHasUncertainty || !isEstCovBlockDiag_)){
	    
	  ::Eigen::MatrixXd *C_inverse_m = node->getInnovationCovInv(); 
	  N = ( C_current - W * *C_inverse_m * W.transpose() ).inverse();	
	  L = (-N * W) * *C_inverse_m;	  
	  K = *C_inverse_m + L.transpose() * N.inverse() * L;
	  double d1 = innov_m->transpose() * L.transpose() * N.inverse() * L * *innov_m;
	  double d2 = 2 * innov.transpose() * L * *innov_m;
	  double d3 = innov.transpose() * N * innov;
	  dist_mahalanobis2 = d2_m + d1 + d2 + d3;

	}else{ // robot pose has no uncertainty or we are dealing with the first measurement

	  dist_mahalanobis2 = d2_m + innov.transpose() * C_current.inverse() * innov;

	}

	// Chi-square gating
	int nDof = measureDim * (1 + nPreviousMeasurements);
	boost::math::chi_squared_distribution<> chi2Dist(nDof);
	double threshold = quantile(chi2Dist, confidence_interval_);
	if(dist_mahalanobis2 < threshold){

	  // gating passed - new node in interpretation tree 
	  JCBB_TreeNodePtr new_node = node->addChild();
	  new_node->setAssociation(z_idx, m_idx, dist_mahalanobis2);
	  int new_innov_size = innov_m->size() + innov.size();
	  new_node->getInnovation()->resize(new_innov_size);
	  *(new_node->getInnovation()) << *innov_m, innov;
	  if(nPreviousMeasurements == 0){

	    *(new_node->getInnovationCovInv()) = C_current.inverse();

	  }else if(isRobotPoseHasUncertainty || !isEstCovBlockDiag_){		
		
	    new_node->getInnovationCovInv()->resize(new_innov_size, new_innov_size);
	    *(new_node->getInnovationCovInv()) << K, L.transpose(), L, N;

	  } 
	    
	  // add to priority queue
	  if(z_idx + 1 < Z_actual_ptr_->size())
	    branchList_.push(new_node);

	  // Keep track of best data association hypothesis
	  if( (nAssociations_best_ < nAssociations + 1) ||  
	      ( nAssociations_best_ == nAssociations + 1 && dist_mahalanobis2_best_ <= dist_mahalanobis2) ){
	    nAssociations_best_ = nAssociations + 1;
	    dist_mahalanobis2_best_ = dist_mahalanobis2;
	    node_best_ = new_node;
	  }
	}

      }

      m_idx++;

    }// BOOST_FOREACH (unused landmark index)
    
    // Add node for outlier measuremnet
    JCBB_TreeNodePtr new_node = node->addChild();
    new_node->setAssociation(z_idx, -1, d2_m);
    *(new_node->getInnovation()) = *(node->getInnovation());
    *(new_node->getInnovationCovInv()) = *(node->getInnovationCovInv());

    if(z_idx + 1 < Z_actual_ptr_->size()){
      branchList_.push(new_node);
    }

    if(  nAssociations_best_ == nAssociations + 1 && dist_mahalanobis2_best_ <= d2_m){
      dist_mahalanobis2_best_ = d2_m;
      node_best_ = new_node;
    }

  }

} // namespace rfs


#endif
