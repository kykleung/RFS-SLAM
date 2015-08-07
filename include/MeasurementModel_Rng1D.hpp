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

#ifndef MEASUREMENTMODEL_RNG1D_HPP
#define MEASUREMENTMODEL_RNG1D_HPP

#include "MeasurementModel.hpp"

namespace rfs{

////////// 1d Measurement Model //////////

/** 
 * \class  MeasurementModel_Rng1D
 * \f[ z = r = h(x, m) + e = m - x + e, e \sim (0, \sigma^2)\f]
 * where 
 * \f$z\f$ or \f$r\f$ is the range measurement, 
 * \f$x\f$ is the robot position,
 * \f$m\f$ is the landmark position
 * \f$e\f$ is the Gaussian noise
 * \brief Range measurement model for 1d point landmarks with Gaussian noise.
 * \author Keith Leung
 */
                                                               
class MeasurementModel_Rng1D : public MeasurementModel <Pose1d, Landmark1d, Measurement1d>{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** \brief Configuration for this MeasurementModel_Rng1D */
  struct Config{
    double probabilityOfDetection_; /**< probability of detection, \f$ P_D \f$ */
    double uniformClutterIntensity_; /**< clutter intensity, \f$ c \f$, assumed to be constant over the sensing range */
    double rangeLimMax_; /**< sensing range limit, beyond which \f$ P_D = 0 \f$*/
    double rangeLimMin_; /**< sensing range limit, below which \f$ P_D = 0 \f$*/
    double rangeLimBuffer_; /**< Used to define a buffer zone around rangeLimMax_ and rangeLimMin_ to indicate a measurement is close to to the sensing limit */
  }config;

 /** Default constructor */
  MeasurementModel_Rng1D();

 /**
  * Constructor that sets the uncertainty of the measurement model, \f$\sigma^2\f$
  * \param[in] Sr range variance \f$\sigma^2\f$
  */
  MeasurementModel_Rng1D(::Eigen::Matrix<double, 1, 1> &Sr);

 /**
  * Constructor that sets the uncertainty of the measurement model, \f$\sigma^2\f$
  * \param[in] Sr range variance \f$\sigma^2\f$
  */
  MeasurementModel_Rng1D(double Sr);

 /** Default destructor */
  ~MeasurementModel_Rng1D();

  /** 
   * Obtain a measurement from a given robot pose and landmark position
   * \f[ z = h(x, m) + e, e \sim (0, \sigma^2) \f] 
   * where \f$z\f$ is a measurement, \f$x\f$ is the robot pose, \f$m\f$ is a landmark position, \f$e\f$ is the zero-mean Gaussian noise.
   * \param[in] pose \f$x\f$, robot pose from which the measurement is made
   * \param[in] landmark \f$m\f$, the measured landmark
   * \param[out] measurement \f$x\f$, the measurement
   * \param[out] jacobian_wrt_lmk if not NULL, the pointed-to matrix is overwritten 
   * by the Jacobian of the measurement model w.r.t. the landmark state evaluated at 
   * \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \param[out] jacobian_wrt_pose if not NULL, the pointed-to matrix is overwritten 
   * by the Jacobian of the measurement model w.r.t. the robot state evaluated at 
   * \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \return true if a valid measurement is produced
   */
  bool measure( const Pose1d &pose, const Landmark1d &landmark, 
		Measurement1d &measurement, 
		::Eigen::Matrix<double, 1, 1> *jacobian_wrt_lmk = NULL,
		::Eigen::Matrix<double, 1, 1> *jacobian_wrt_pose = NULL);

  /** 
   * \f[ m = h^{-1}(x, z)\f] 
   * where \f$z\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
   * \param[in] pose \f$x\f$, robot pose (the uncertainty is not used here because the 
   * RBPHDFilter represents robot pose estimates with particles)
   * \param[in] measurement  \f$z\f$ measurement, for which the uncertainty is \f$\sigma^2\f$
   * \param[out] landmark  \f$m\f$, predicted landmark position with uncertainty
   */
  void inverseMeasure(const Pose1d &pose, const Measurement1d &measurement, Landmark1d &landmark);

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
  double probabilityOfDetection( const Pose1d &pose,
				 const Landmark1d &landmark,
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

}

#endif
