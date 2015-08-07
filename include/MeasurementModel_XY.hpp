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


#ifndef MEASUREMENTMODEL_XY_HPP
#define MEASUREMENTMODEL_XY_HPP

#include "MeasurementModel.hpp"

namespace rfs{

////////// 2d Range-Bearing Measurement Model //////////

/** 
 * \class MeasurementModel_XY
 * A range and bearing measurement model for 2d point landmarks, with Gaussian noise.
 * \f[ \mathbf{z} = 
 * \begin{bmatrix} z_x \\ z_y \end{bmatrix} =
 * \mathbf{h}(\mathbf{x}, \mathbf{m}) + \mathbf{e} = 
 * \mathbf{h}\left(\begin{bmatrix}x \\ y \\ \theta\end{bmatrix}, \begin{bmatrix}x_m \\ y_m\end{bmatrix}\right) + \mathbf{e} = 
 * \begin{bmatrix} \cos\theta (x_m - x) + \sin\theta (y_m - y) \\ -\sin\theta (x_m - x) + \cos\theta (y_m - y) \end{bmatrix} + \mathbf{e} , \quad \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f]
 * where
 * \f$\mathbf{z} = (z_x, z_y)\f$ is the x-y measurement relative to the robot,
 * \f$\mathbf{x} = (x, y, \theta)\f$ is the robot pose,
 * \f$\mathbf{m} = (x_m, y_m)\f$ is the landmark position, and
 * \f$\mathbf{e}\f$ is the noise with covariance \f$\mathbf{R}\f$
 * \brief A 2d x-y measurement model
 * \author Felipe Inostroza, Keith Leung 
 */
                                                               
class MeasurementModel_XY : public MeasurementModel <Pose2d, Landmark2d, Measurement2d>{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** \brief Configuration for this 2d MeasurementModel_XY */
  struct Config{
    double probabilityOfDetection_; /**< probability of detection, \f$ P_D \f$ */
    double uniformClutterIntensity_; /**< clutter intensity, \f$ c \f$, assumed to be constant over the sensing area */
    double rangeLimMax_; /**< sensing range limit, beyond which \f$ P_D = 0 \f$*/
    double rangeLimMin_; /**< sensing range limit, below which \f$ P_D = 0 \f$*/
    double rangeLimBuffer_; /**< Used to define a buffer zone around rangeLimMax_ and rangeLimMin_ to indicate a measurement is close to to the sensing limit */
  }config;

 /** Default constructor */
  MeasurementModel_XY();

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model, \f$\mathbf{R}\f$
  * \param covZ measurement covariance, \f$\mathbf{R}\f$
  */
  MeasurementModel_XY(::Eigen::Matrix2d &covZ);

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model, 
  * \f[\mathbf{R} = \begin{bmatrix} \sigma_{z_x}^2 & 0 \\ 0 & \sigma_{z_y}^2 \end{bmatrix}\f]
  * range and bearing are assumed to be uncorrelated
  * \param Sx x variance \f$\sigma_{z_x}^2\f$
  * \param Sy y variance \f$\sigma_{z_y}^2\f$
  */
  MeasurementModel_XY(double Sx, double Sy);

 /** Default destructor */
  ~MeasurementModel_XY();

  /** 
   * Obtain a measurement from a given robot pose and landmark position
   * \f[ \mathbf{z} = \mathbf{h}(\mathbf{x}, \mathbf{m} ) + \mathbf{e}, \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f] 
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position, \f$\mathbf{e}\f$ is the zero-mean Gaussian noise.
   * \param[in] pose \f$\mathbf{x}\f$, robot pose from which the measurement is made
   * \param[in] landmark \f$\mathbf{m}\f$, the measured landmark
   * \param[out] measurement \f$\mathbf{x}\f$, the measurement
   * \param[out] jacobian_wrt_lmk if not NULL, the pointed-to matrix is overwritten 
   * by the Jacobian of the measurement model w.r.t. the landmark state evaluated at 
   * \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \param[out] jacobian_wrt_pose if not NULL, the pointed-to matrix is overwritten 
   * by the Jacobian of the measurement model w.r.t. the robot state evaluated at 
   * \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$.
   * \return true if a valid measurement is produced
   */
  bool measure( const Pose2d &pose, const Landmark2d &landmark, 
		Measurement2d &measurement, 
		::Eigen::Matrix2d *jacobian_wrt_lmk = NULL,
		::Eigen::Matrix<double, 2, 3> *jacobian_wrt_pose = NULL);

  /** 
   * \f[ \mathbf{m} = \mathbf{h}^{-1}(\mathbf{x}, \mathbf{z} )\f] 
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
   * \param[in] pose \f$\mathbf{x}\f$, robot pose
   * \param[in] measurement  \f$\mathbf{z}\f$ measurement, for which the uncertainty is \f$\mathbf{R}\f$
   * \param[out] landmark  \f$\mathbf{m}\f$, predicted landmark position with uncertainty
   */
  void inverseMeasure(const Pose2d &pose, const Measurement2d &measurement, Landmark2d &landmark);

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
  double probabilityOfDetection( const Pose2d &pose,
				 const Landmark2d &landmark,
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


}

#endif
