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

#include "MeasurementModel_Rng1D.hpp"

namespace rfs
{

MeasurementModel_Rng1D::MeasurementModel_Rng1D(){

  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_Rng1D::MeasurementModel_Rng1D(Eigen::Matrix<double, 1, 1> &Sr){
  setNoise(Sr);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_Rng1D::MeasurementModel_Rng1D(double Sr){
  Measurement1d::Mat S;
  S << Sr;
  setNoise(S);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_Rng1D::~MeasurementModel_Rng1D(){}

bool MeasurementModel_Rng1D::measure(const Pose1d &pose, 
				 const Landmark1d &landmark, 
				 Measurement1d &measurement, 
				 Eigen::Matrix<double, 1, 1> *jacobian){

  Pose1d::Vec robotPos;
  Landmark1d::Vec lmkPos;
  Landmark1d::Mat lmkPosUncertainty;
  Measurement1d::Vec z;
  Eigen::Matrix<double, 1, 1> H;
  Eigen::Matrix<double, 1, 1> var;

  pose.get(robotPos);
  landmark.get(lmkPos, lmkPosUncertainty);
  z = lmkPos - robotPos;
  H << 1;
  var = H * lmkPosUncertainty * H.transpose() + R_;

  measurement.set(z, var);

  if(jacobian != NULL)
    *jacobian = H;

  if(fabs(z(0)) > config.rangeLimMax_ || fabs(z(0)) < config.rangeLimMin_)
    return false;
  else
    return true;

}

void MeasurementModel_Rng1D::inverseMeasure(const Pose1d &pose, 
					const Measurement1d &measurement, 
					Landmark1d &landmark){

  Pose1d::Vec x;
  Measurement1d::Vec z;
  Landmark1d::Vec m;
 
  pose.get(x);
  measurement.get(z);
  landmark.get(m);

  m = x + z;
  landmark.set( m, R_ );

}

double MeasurementModel_Rng1D::probabilityOfDetection( const Pose1d &pose,
						   const Landmark1d &landmark,
						   bool &isCloseToSensingLimit ){

  Pose1d::Vec robotPose;
  Landmark1d::Vec landmarkState;
  double range, Pd;
  
  isCloseToSensingLimit = false;

  pose.get(robotPose);
  landmark.get(landmarkState);
  range = fabs( landmarkState(0) - robotPose(0) );

  if( range <= config.rangeLimMax_ && range >= config.rangeLimMin_){
    Pd = config.probabilityOfDetection_;
  }else{
    Pd = 0;
  } 

  if( ( range >= (config.rangeLimMax_ - config.rangeLimBuffer_ ) &&
	range <= (config.rangeLimMax_ + config.rangeLimBuffer_ ) ) ||
      ( range >= (config.rangeLimMin_ - config.rangeLimBuffer_ ) &&
	range <= (config.rangeLimMin_ + config.rangeLimBuffer_ ) ) )
    isCloseToSensingLimit = true;

  return Pd;
}

double MeasurementModel_Rng1D::clutterIntensity( Measurement1d &z,
					     int nZ){
  return config.uniformClutterIntensity_;
}


double MeasurementModel_Rng1D::clutterIntensityIntegral( int nZ ){
  double sensingLength = config.rangeLimMax_ - config.rangeLimMin_;
  return ( config.uniformClutterIntensity_ * sensingLength );
}

}
