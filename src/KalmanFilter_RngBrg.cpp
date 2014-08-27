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

#include "KalmanFilter_RngBrg.hpp"

namespace rfs
{

KalmanFilter_RngBrg::KalmanFilter_RngBrg(){
  config.rangeInnovationThreshold_ = -1;
  config.bearingInnovationThreshold_ = -1;
};
  
KalmanFilter_RngBrg::KalmanFilter_RngBrg(StaticProcessModel<Landmark2d> *pProcessModel,
						   MeasurementModel_RngBrg *pMeasurementModel):
  KalmanFilter<StaticProcessModel<Landmark2d>, MeasurementModel_RngBrg>
  (pProcessModel, pMeasurementModel){
  config.rangeInnovationThreshold_ = -1;
  config.bearingInnovationThreshold_ = -1;
}

bool KalmanFilter_RngBrg::calculateInnovation(Vec &z_exp, Vec &z_act, Vec &z_innov){
    
  z_innov = z_act - z_exp;
  
  if(config.rangeInnovationThreshold_ > 0 && fabs(z_innov(0)) > config.rangeInnovationThreshold_)
    return false;
  while(z_innov(1) > PI)
    z_innov(1) -= 2 * PI;
  while(z_innov(1) < -PI)
    z_innov(1) += 2 * PI;
  if(config.bearingInnovationThreshold_ > 0 && fabs(z_innov(1) > config.bearingInnovationThreshold_) )
    return false;
  return true;
}
  
}
