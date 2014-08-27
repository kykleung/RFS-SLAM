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

// Class for using the Kalman filter equations
// Felipe Inostroza 2013, Keith Leung 2014

#ifndef KALMANFILTER_RNGBRG_HPP
#define KALMANFILTER_RNGBRG_HPP

#include "KalmanFilter.hpp"
#include "MeasurementModel_RngBrg.hpp"

namespace rfs{

/**
 * \class KalmanFilter_RngBrg
 * \brief A Kalman filter for updating a 2d landmark position from 
 * a single range-bearing measurements. This is derived from the base
 * Kalman Filter to handle innovation involving a rotation (bearing)
 */
class KalmanFilter_RngBrg : public KalmanFilter <StaticProcessModel<Landmark2d>, MeasurementModel_RngBrg>{

  typedef MeasurementModel_RngBrg::TMeasurement::Vec Vec;

public:

  /** \brief Configuration for this KalmanFilter_RngBrg */
  struct Config{
    /** If positive, the innovation threshold above which an update is not processed for stability reasons. */
    double rangeInnovationThreshold_; 
    /** If positive, the innovation threshold above which an update is not processed for stability reasons. */
    double bearingInnovationThreshold_;
  }config;

  /**
   * Default constructor
   */
  KalmanFilter_RngBrg();
  
  /**
   * Constructor
   * \param[in] pMeasurementModel Pointer to the measurement model
   * \param[in] pProcessModel Pointer to the process model
   */
  KalmanFilter_RngBrg(StaticProcessModel<Landmark2d> *pProcessModel,
			   MeasurementModel_RngBrg *pMeasurementModel);

  /**
   * Function to calculate the innovation 
   * \param[in] z_exp expected measurement predicted using the measurement model 
   * \param[in] z_act actual measurement
   * \param[out] z_innov innovation
   */
  bool calculateInnovation(Vec &z_exp, Vec &z_act, Vec &z_innov);

};

}
#endif
