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

// Measurement model
// Felipe Inostroza, Keith Leung 2013

#ifndef MEASUREMENTMODEL_HPP
#define MEASUREMENTMODEL_HPP

#include "Measurement.hpp"
#include "Landmark.hpp"
#include "Pose.hpp"

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
  void getNoise( typename MeasurementType::Mat &R ){
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
   * \param[out] jacobian if not NULL, the pointed-to matrix is overwritten 
   * by the Jacobian of the measurement model, \f$\mathbf{H}\f$, evaluated at \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \return true if a valid measurement is produced
   */
  virtual bool measure( PoseType &pose, LandmarkType &landmark, 
			MeasurementType &meaurement, 
			Eigen::Matrix<double , 
				      MeasurementType::Vec::RowsAtCompileTime ,
				      LandmarkType::Vec::RowsAtCompileTime > 
			*jacobian = NULL ) = 0;
  
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
	
    typename MeasurementType::Vec z;
    typename MeasurementType::Vec noise;

    PoseType* pose_sample = &pose;
    LandmarkType* landmark_sample = &landmark;
    bool deallocatePose = false;
    bool deallocateLandmark = false;
    
    if(usePoseWhiteGaussianNoise){

      pose_sample = new PoseType;
      deallocatePose = true;      
      pose.sample(*pose_sample);
    }

    if(useLandmarkWhiteGaussianNoise){

      landmark_sample = new LandmarkType;
      deallocateLandmark = true;
      landmark.sample(*landmark_sample);
    }

    bool success = this->measure( *pose_sample, *landmark_sample, measurement);

    if(success){
      if(useAdditiveWhiteGaussianNoise){
	measurement.setCov(R_);
	measurement.sample();
      }
    }

    if(deallocatePose)
      delete pose_sample;
    if(deallocateLandmark)
      delete landmark_sample;

    return success;

  }

  /** 
   * Abstract function for the inverse measurement model
   * \f[ \mathbf{m} = \mathbf{h}^{-1}(\mathbf{x}, \mathbf{z} )\f] 
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
   * \note This must be implemented in a derived class, and both the mean and the covariance of \f$\mathbf{m}\f$ should 
   * be calculated. This is used in the RBPHDFilter for generating birth Gaussians
   * \param[in] pose \f$\mathbf{x}\f$, robot pose (the uncertainty is not used here because the 
   * RBPHDFilter represents robot pose estimates with particles)
   * \param[in] measurement  \f$\mathbf{z}\f$ measurement, for which the uncertainty is \f$\mathbf{R}\f$
   * \param[out] landmark  \f$\mathbf{m}\f$, predicted landmark position with uncertainty
   */
  virtual void inverseMeasure( PoseType &pose,
			       MeasurementType &measurement, 
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
  virtual double probabilityOfDetection( PoseType &pose,
					 LandmarkType &landmark,
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


////////// 2d Range-Bearing Measurement Model //////////

/** 
 * \class RangeBearingModel
 * A range and bearing measurement model for 2d point landmarks, with Gaussian noise.
 * \f[ \mathbf{z} = 
 * \begin{bmatrix} r \\ b \end{bmatrix} =
 * \mathbf{h}(\mathbf{x}, \mathbf{m}) + \mathbf{e} = 
 * \mathbf{h}\left(\begin{bmatrix}x \\ y \\ \theta\end{bmatrix}, \begin{bmatrix}x_m \\ y_m\end{bmatrix}\right) + \mathbf{e} = 
 * \begin{bmatrix} \sqrt{(x_m - x)^2+(y_m - y)^2}) \\ \arctan{\left(\frac{y_m - y}{x_m - x}\right)} - \theta \end{bmatrix} + \mathbf{e} , \quad \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f]
 * where
 * \f$\mathbf{z} = (r, b)\f$ is the range and bearing measurement,
 * \f$\mathbf{x} = (x, y, \theta)\f$ is the robot pose,
 * \f$\mathbf{m} = (x_m, y_m)\f$ is the landmark position, and
 * \f$\mathbf{e}\f$ is the noise with covariance \f$\mathbf{R}\f$
 * \brief A 2d range-bearing measurement model
 * \author Felipe Inostroza, Keith Leung 
 */
                                                               
class RangeBearingModel : public MeasurementModel <Pose2d, Landmark2d, Measurement2d>{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** \brief Configuration for this 2d RangeBearingModel */
  struct Config{
    double probabilityOfDetection_; /**< probability of detection, \f$ P_D \f$ */
    double uniformClutterIntensity_; /**< clutter intensity, \f$ c \f$, assumed to be constant over the sensing area */
    double rangeLimMax_; /**< sensing range limit, beyond which \f$ P_D = 0 \f$*/
    double rangeLimMin_; /**< sensing range limit, below which \f$ P_D = 0 \f$*/
    double rangeLimBuffer_; /**< Used to define a buffer zone around rangeLimMax_ and rangeLimMin_ to indicate a measurement is close to to the sensing limit */
  }config;

 /** Default constructor */
  RangeBearingModel();

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model, \f$\mathbf{R}\f$
  * \param covZ measurement covariance, \f$\mathbf{R}\f$
  */
  RangeBearingModel(Eigen::Matrix2d &covZ);

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model, 
  * \f[\mathbf{R} = \begin{bmatrix} \sigma_r^2 & 0 \\ 0 & \sigma_b^2 \end{bmatrix}\f]
  * range and bearing are assumed to be uncorrelated
  * \param Sr Range variance \f$\sigma_r^2\f$
  * \param Sb Bearing variance \f$\sigma_b^2\f$
  */
  RangeBearingModel(double Sr, double Sb);

 /** Default destructor */
  ~RangeBearingModel();

  /** 
   * Obtain a measurement from a given robot pose and landmark position
   * \f[ \mathbf{z} = \mathbf{h}(\mathbf{x}, \mathbf{m} ) + \mathbf{e}, \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f] 
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position, \f$\mathbf{e}\f$ is the zero-mean Gaussian noise.
   * \param[in] pose \f$\mathbf{x}\f$, robot pose from which the measurement is made
   * \param[in] landmark \f$\mathbf{m}\f$, the measured landmark
   * \param[out] measurement \f$\mathbf{x}\f$, the measurement
   * \param[out] jacobian if not NULL, the pointed-to matrix is overwritten 
   * by the Jacobian of the measurement model, \f$\mathbf{H}\f$, evaluated at \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \return true if a valid measurement is produced
   */
  bool measure( Pose2d &pose, Landmark2d &landmark, 
		Measurement2d &measurement, Eigen::Matrix2d *jacobian = NULL);

  /** 
   * \f[ \mathbf{m} = \mathbf{h}^{-1}(\mathbf{x}, \mathbf{z} )\f] 
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
   * \param[in] pose \f$\mathbf{x}\f$, robot pose (the uncertainty is not used here because the 
   * RBPHDFilter represents robot pose estimates with particles)
   * \param[in] measurement  \f$\mathbf{z}\f$ measurement, for which the uncertainty is \f$\mathbf{R}\f$
   * \param[out] landmark  \f$\mathbf{m}\f$, predicted landmark position with uncertainty
   */
  void inverseMeasure(Pose2d &pose, Measurement2d &measurement, Landmark2d &landmark);

  /**
   * Abstract function of determining a landmark's probability of detection, and if the landmark is close to the sensing limit.
   * Through this we can indirectly specify sensing limits and other sensor characteristics
   * The probability of detection is necessary as a parameter is the PHD Filter. Indicating whether a landmark is close to the 
   * sensing limit matters in the implementation for providing a better map estimate, as it reduces landmark disappearance
   * near the sensing limit due to the probability of detection mismatch.
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] isCloseToSensingLimit true if landmark is close to the sensing limit
   * \return probability of detection
   */
  double probabilityOfDetection( Pose2d &pose,
				 Landmark2d &landmark,
				 bool &isCloseToSensingLimit);

  /**
   * Determine the clutter intensity in measurement space.
   * Uniform clutter intensity is assumed
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity
   */
  double clutterIntensity( Measurement2d &z,
			   int nZ );

  /**
   * Determine the clutter intensity integral in measurement space.
   * This is calculated based on the probablity of false alarm,
   * defined as p( NULL | measurement exists)
   * \param[in] nZ the cardinality of Z
   * \return clutter intensity
   */
  double clutterIntensityIntegral( int nZ = 0);

};



////////// 1d Measurement Model //////////

/** 
 * \class  MeasurementModel1d
 * \f[ z = r = h(x, m) + e = m - x + e, e \sim (0, \sigma^2)\f]
 * where 
 * \f$z\f$ or \f$r\f$ is the range measurement, 
 * \f$x\f$ is the robot position,
 * \f$m\f$ is the landmark position
 * \f$e\f$ is the Gaussian noise
 * \brief Range measurement model for 1d point landmarks with Gaussian noise.
 * \author Keith Leung
 */
                                                               
class MeasurementModel1d : public MeasurementModel <Pose1d, Landmark1d, Measurement1d>{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** \brief Configuration for this MeasurementModel1d */
  struct Config{
    double probabilityOfDetection_; /**< probability of detection, \f$ P_D \f$ */
    double uniformClutterIntensity_; /**< clutter intensity, \f$ c \f$, assumed to be constant over the sensing range */
    double rangeLimMax_; /**< sensing range limit, beyond which \f$ P_D = 0 \f$*/
    double rangeLimMin_; /**< sensing range limit, below which \f$ P_D = 0 \f$*/
    double rangeLimBuffer_; /**< Used to define a buffer zone around rangeLimMax_ and rangeLimMin_ to indicate a measurement is close to to the sensing limit */
  }config;

 /** Default constructor */
  MeasurementModel1d();

 /**
  * Constructor that sets the uncertainty of the measurement model, \f$\sigma^2\f$
  * \param[in] Sr range variance \f$\sigma^2\f$
  */
  MeasurementModel1d(Eigen::Matrix<double, 1, 1> &Sr);

 /**
  * Constructor that sets the uncertainty of the measurement model, \f$\sigma^2\f$
  * \param[in] Sr range variance \f$\sigma^2\f$
  */
  MeasurementModel1d(double Sr);

 /** Default destructor */
  ~MeasurementModel1d();

  /** 
   * Obtain a measurement from a given robot pose and landmark position
   * \f[ z = h(x, m) + e, e \sim (0, \sigma^2) \f] 
   * where \f$z\f$ is a measurement, \f$x\f$ is the robot pose, \f$m\f$ is a landmark position, \f$e\f$ is the zero-mean Gaussian noise.
   * \param[in] pose \f$x\f$, robot pose from which the measurement is made
   * \param[in] landmark \f$m\f$, the measured landmark
   * \param[out] measurement \f$x\f$, the measurement
   * \param[out] jacobian if not NULL, the pointed-to matrix is overwritten 
   * by the derivative of the measurement model, \f$\mathbf{H}\f$, evaluated at \f$x\f$ and \f$m\f$
   * \return true if a valid measurement is produced
   */
  bool measure( Pose1d &pose, Landmark1d &landmark, 
		Measurement1d &measurement, Eigen::Matrix<double, 1, 1> *jacobian = NULL);

  /** 
   * \f[ m = h^{-1}(x, z)\f] 
   * where \f$z\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
   * \param[in] pose \f$x\f$, robot pose (the uncertainty is not used here because the 
   * RBPHDFilter represents robot pose estimates with particles)
   * \param[in] measurement  \f$z\f$ measurement, for which the uncertainty is \f$\sigma^2\f$
   * \param[out] landmark  \f$m\f$, predicted landmark position with uncertainty
   */
  void inverseMeasure(Pose1d &pose, Measurement1d &measurement, Landmark1d &landmark);

  /**
   * Abstract function of determining a landmark's probability of detection, and if the landmark is close to the sensing limit.
   * Through this we can indirectly specify sensing limits and other sensor characteristics
   * The probability of detection is necessary as a parameter is the PHD Filter. Indicating whether a landmark is close to the 
   * sensing limit matters in the implementation for providing a better map estimate, as it reduces landmark disappearance
   * near the sensing limit due to the probability of detection mismatch.
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] isCloseToSensingLimit true if landmark is close to the sensing limit
   * \return probability of detection
   */
  double probabilityOfDetection( Pose1d &pose,
				 Landmark1d &landmark,
				 bool &isCloseToSensingLimit);

  /**
   * Determine the clutter intensity in measurement space.
   * Uniform clutter intensity is assumed
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity
   */
  double clutterIntensity( Measurement1d &z,
			   int nZ );

  /**
   * Determine the clutter intensity integral in measurement space.
   * This is calculated based on the probablity of false alarm,
   * defined as p( NULL | measurement exists)
   * \param[in] nZ the cardinality of Z
   * \return clutter intensity
   */
  double clutterIntensityIntegral( int nZ = 0);

};



#endif
