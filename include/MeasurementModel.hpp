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

#ifndef MEASUREMENTMODEL_HPP
#define MEASUREMENTMODEL_HPP

#include "Measurement.hpp"
#include "Landmark.hpp"
#include "Pose.hpp"

namespace rfs{

/** 
 * \class MeasurementModel
 * \brief An abstract class for defining the measurement model
 *
 * \f[ \mathbf{z} = \mathbf{h}(\mathbf{x}, \mathbf{m} ) + \mathbf{e}, \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f] 
 * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position, \f$\mathbf{e}\f$ is the zero-mean Gaussian noise.
 * \tparam PoseType sensor pose type
 * \tparam LandmarkType measured object type
 * \tparam MeasurementType measurement type
 * \author Felipe Inostroza, Keith Leung
 */
template<class PoseType, class LandmarkType, class MeasurementType>
class MeasurementModel
{
public:

  typedef PoseType TPose;
  typedef LandmarkType TLandmark;
  typedef MeasurementType TMeasurement;
  typedef ::Eigen::Matrix<double, 
			  MeasurementType::Vec::RowsAtCompileTime ,
			  LandmarkType::Vec::RowsAtCompileTime> TJacobianLmk;
  typedef ::Eigen::Matrix<double,
			  MeasurementType::Vec::RowsAtCompileTime,
			  PoseType::Vec::RowsAtCompileTime> TJacobianPose;
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** Default constructor */
  MeasurementModel() : R_( MeasurementType::Mat::Zero()) {}

  /** Default destructor */
  ~MeasurementModel(){}

  /** 
   * Set the zero-mean-white-Gaussian additive noise covariance matrix, \f$\mathbf{R}\f$
   * \param[in] R covariance matrix
   */
  void setNoise( typename MeasurementType::Mat &R ){
    R_ = R;
  }

  /** 
   * Get the zero-mean-white-Gaussian noise covariance matrix, \f$\mathbf{R}\f$
   * \param[out] R covariance matrix
   */
  void getNoise( typename MeasurementType::Mat &R ) const{
    R = R_;
  }
  
  /** 
   * Abstract function for predicting a measurement from a robot pose and a landmark position
   * \f[ \mathbf{z} = \mathbf{h}(\mathbf{x}, \mathbf{m} ) + \mathbf{e}, \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f] 
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position, \f$\mathbf{e}\f$ is the zero-mean Gaussian noise.
   * \note This must be implemented in a derived class
   * \param[in] pose \f$\mathbf{x}\f$, robot pose from which the measurement is made
   * \param[in] landmark \f$\mathbf{m}\f$, the measured landmark
   * \param[out] measurement \f$\mathbf{x}\f$, the measurement
   * \param[out] jacobian_wrt_lmk if not NULL, the pointed-to matrix is overwritten 
   * by the Jacobian of the measurement model w.r.t. the landmark state evaluated at 
   * \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \param[out] jacobian_wrt_pose if not NULL, the pointed-to matrix is overwritten 
   * by the Jacobian of the measurement model w.r.t. the robot state evaluated at 
   * \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \return true if a valid measurement is produced
   */
  virtual bool measure( const PoseType &pose, 
			const LandmarkType &landmark, 
			MeasurementType &measurement, 
			::Eigen::Matrix<double, 
					MeasurementType::Vec::RowsAtCompileTime ,
					LandmarkType::Vec::RowsAtCompileTime> *jacobian_wrt_lmk = NULL,
			::Eigen::Matrix<double,
			                MeasurementType::Vec::RowsAtCompileTime,
					PoseType::Vec::RowsAtCompileTime> *jacobian_wrt_pose = NULL 
			) = 0;
  
  /**
   * Sample a measurement with noise parameters from the model, robot pose, and landmark position.
   * \param[in] pose \f$\mathbf{x}\f$, robot pose from which the measurement is made
   * \param[in] landmark \f$\mathbf{m}\f$, the measured landmark
   * \param[out] measurement sampled measurement (which does not contain uncertainty information, i.e., the covariance information is not set for this RandomVec)
   * \param[in] useAdditiveWhiteGaussianNoise if true, sample from the model's noise covariance, \f$\mathbf{R}\f$
   * \param[in] usePoseWhiteGaussianNoise if true, sample the uncertain robot pose using its covariance, 
   * \f$\boldsymbol{\Sigma}_{\mathbf{x}}\f$, and interpret as zero-mean-white-Gaussian noise
   * \param[in] useLandmarkWhiteGaussianNoise if true, sample the uncertain landmark position using its covariance,  
   * \f$\boldsymbol{\Sigma}_{\mathbf{m}}\f$, and interpret as zero-mean-white-Gaussian noise
   * \return true if sucessfully sampled a measurement
   */
  bool sample( PoseType &pose, LandmarkType &landmark, 
	       MeasurementType &measurement,
      	       bool useAdditiveWhiteGaussianNoise = true,		       
	       bool usePoseWhiteGaussianNoise = false,
	       bool useLandmarkWhiteGaussianNoise = false){
	      
    PoseType pose_sample; /**< Sampled pose*/
    LandmarkType landmark_sample; /**< Sampled landmark*/
    
    if(usePoseWhiteGaussianNoise){
      pose.sample(pose_sample);
    }else{
      pose_sample = pose;
    }

    if(useLandmarkWhiteGaussianNoise){
      landmark.sample(landmark_sample);
    }else{
      landmark_sample = landmark;
    }

    bool success = this->measure( pose_sample, landmark_sample, measurement);

    if(success && useAdditiveWhiteGaussianNoise){
      measurement.setCov(R_);
      measurement.sample();
    }

    return success;

  }

  /** 
   * Abstract function for the inverse measurement model
   * \f[ \mathbf{m} = \mathbf{h}^{-1}(\mathbf{x}, \mathbf{z} )\f] 
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
   * \note This must be implemented in a derived class, and both the mean and the covariance of \f$\mathbf{m}\f$ should 
   * be calculated. This is used in the RBPHDFilter for generating birth Gaussians
   * \param[in] pose \f$\mathbf{x}\f$, robot pose (the uncertainty is not used here because the pose_sample_
   * RBPHDFilter represents robot pose estimates with particles)
   * \param[in] measurement  \f$\mathbf{z}\f$ measurement, for which the uncertainty is \f$\mathbf{R}\f$
   * \param[out] landmark  \f$\mathbf{m}\f$, predicted landmark position with uncertainty
   */
  virtual void inverseMeasure( const PoseType &pose,
			       const MeasurementType &measurement, 
			       LandmarkType &landmark ) = 0;

  /**
   * Abstract function of determining a landmark's probability of detection, and if the landmark is close to the sensing limit.
   * Through this we can indirectly specify sensing limits and other sensor characteristics
   * The probability of detection is necessary as a parameter is the PHD Filter. Indicating whether a landmark is close to the 
   * sensing limit matters in the implementation for providing a better map estimate, as it reduces landmark disappearance
   * near the sensing limit due to the probability of detection mismatch.
   * \note If this is not reimplemented in a derived class, it will always
   * return a probability of detection of 1
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] isCloseToSensingLimit true if landmark is close to the sensing limit
   * \return probability of detection
   */
  virtual double probabilityOfDetection( const PoseType &pose,
					 const LandmarkType &landmark,
					 bool &isCloseToSensingLimit){
    isCloseToSensingLimit = false;
    return 1;
  }

  /**
   * Abstract function for determining the clutter intensity, \f$c\f$
   * \note This should be reimplemented in a derived class
   * \param[in] z measurement point at which clutter intensity is to be determined
   * \param[in] nZ the cardinality of measurement set Z, of which z is a member.
   * \return clutter intensity
   */
  virtual double clutterIntensity( MeasurementType &z,
				   int nZ ){
    return 0;
  }

  /**
   * Abstract function for determining the clutter intensity integral over the measurement space
   * (i.e., the expected number of clutter measurements) 
   * \note This should be reimplemented in a derived class
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity integral
   */
  virtual double clutterIntensityIntegral( int nZ ){
    return 0;
  }  

protected:

  typename MeasurementType::Mat R_; /**< additive zero-mean Gaussian noise covariance */





};

}
#endif
