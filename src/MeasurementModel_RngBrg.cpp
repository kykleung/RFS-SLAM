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

#include "MeasurementModel_RngBrg.hpp"

namespace rfs
{

MeasurementModel_RngBrg::MeasurementModel_RngBrg(){

  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}


MeasurementModel_RngBrg::MeasurementModel_RngBrg(Eigen::Matrix2d &covZ){

  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_RngBrg::MeasurementModel_RngBrg(double Sr, double Sb){

  Eigen::Matrix2d covZ;
  covZ <<  Sr, 0, 0, Sb;
  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_RngBrg::~MeasurementModel_RngBrg(){}

bool MeasurementModel_RngBrg::measure(const Pose2d &pose, 
				const Landmark2d &landmark, 
				Measurement2d &measurement, 
				Eigen::Matrix2d *jacobian){

  Eigen::Vector3d robotPose;
  Eigen::Vector2d mean, landmarkState;
  Eigen::Matrix2d H, landmarkUncertainty, cov;
  double range, bearing;

  pose.get(robotPose);
  landmark.get(landmarkState,landmarkUncertainty);

  range = sqrt(  pow(landmarkState(0) - robotPose(0), 2)
		+pow(landmarkState(1) - robotPose(1), 2) );

  bearing = atan2( landmarkState(1) - robotPose(1) , landmarkState(0) - robotPose(0) ) - robotPose(2);

  while(bearing>PI) bearing-=2*PI;
  while(bearing<-PI) bearing+=2*PI;

  mean << range, bearing ;

  H <<  (landmarkState(0)-robotPose(0))/mean(0)          , (landmarkState(1)-robotPose(1))/mean(0) ,
        -(landmarkState(1)-robotPose(1))/(pow(mean(0),2)), (landmarkState(0)-robotPose(0))/pow(mean(0),2) ;
  
  cov = H * landmarkUncertainty * H.transpose() + R_;
  measurement.set(mean, cov);

  
  if(jacobian != NULL)
    *jacobian = H;

  if(range > config.rangeLimMax_ || range < config.rangeLimMin_)
    return false;
  else
    return true;
}

void MeasurementModel_RngBrg::inverseMeasure(const Pose2d &pose, 
				       const Measurement2d &measurement, 
				       Landmark2d &landmark){
  Eigen::Vector3d poseState;
  Eigen::Vector2d measurementState, mean;
  Eigen::Matrix2d measurementUncertainty, covariance, Hinv;
 
  pose.get(poseState);
  measurement.get(measurementState);
  this->getNoise(measurementUncertainty); 
  mean << poseState(0) + measurementState(0) *cos( poseState(2) + measurementState(1) ),
          poseState(1) + measurementState(0) *sin( poseState(2) + measurementState(1) );

  Hinv << cos(poseState(2)+measurementState(1)) , -measurementState(0)*sin(poseState(2)+measurementState(1)) ,
          sin(poseState(2)+measurementState(1)) , measurementState(0)*cos(poseState(2)+measurementState(1));

  covariance = Hinv * measurementUncertainty *Hinv.transpose();
  landmark.set( mean, covariance );

}

double MeasurementModel_RngBrg::probabilityOfDetection( const Pose2d &pose,
						  const Landmark2d &landmark,
						  bool &isCloseToSensingLimit ){

  Pose2d::Vec robotPose;
  Landmark2d::Vec landmarkState;
  double range, Pd;
  
  isCloseToSensingLimit = false;

  pose.get(robotPose);
  landmark.get(landmarkState);

  range = sqrt(  pow(landmarkState(0) - robotPose(0), 2)
		+pow(landmarkState(1) - robotPose(1), 2) );

  if( range <= config.rangeLimMax_ && range >= config.rangeLimMin_){
    Pd = config.probabilityOfDetection_;
    if( range >= (config.rangeLimMax_ - config.rangeLimBuffer_ ) ||
	range <= (config.rangeLimMin_ + config.rangeLimBuffer_ ) )
      isCloseToSensingLimit = true;
  }else{
    Pd = 0;
    if( range <= (config.rangeLimMax_ + config.rangeLimBuffer_ ) &&
	range >= (config.rangeLimMin_ - config.rangeLimBuffer_ ) )
      isCloseToSensingLimit = true;
  } 

  return Pd;
}

double MeasurementModel_RngBrg::clutterIntensity( Measurement2d &z,
					    int nZ){
  return config.uniformClutterIntensity_;
}


double MeasurementModel_RngBrg::clutterIntensityIntegral( int nZ ){
  double sensingArea_ = 2 * PI * (config.rangeLimMax_ - config.rangeLimMin_);
  return ( config.uniformClutterIntensity_ * sensingArea_ );
}

}
