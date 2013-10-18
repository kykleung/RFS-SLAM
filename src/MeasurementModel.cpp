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

#include <Eigen/LU>
#include <math.h>
#include "MeasurementModel.hpp"

/******** Implementation of 2d measurement Model (Range and Bearing) **********/

RangeBearingModel::RangeBearingModel(){

  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}


RangeBearingModel::RangeBearingModel(Eigen::Matrix2d &covZ){

  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

RangeBearingModel::RangeBearingModel(double Sr, double Sb){

  Eigen::Matrix2d covZ;
  covZ <<  Sr, 0, 0, Sb;
  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

RangeBearingModel::~RangeBearingModel(){}

bool RangeBearingModel::measure(Pose2d  &pose, Landmark2d &landmark, 
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

void RangeBearingModel::inverseMeasure(Pose2d &pose, Measurement2d &measurement, 
				       Landmark2d &landmark){
  Eigen::Vector3d poseState;
  Eigen::Vector2d measurementState, mean;
  Eigen::Matrix2d measurementUncertainty, covariance, Hinv;
  double t;
 
  pose.get(poseState);
  measurement.get(measurementState, t);
  this->getNoise(measurementUncertainty); 
  mean << poseState(0) + measurementState(0) *cos( poseState(2) + measurementState(1) ),
          poseState(1) + measurementState(0) *sin( poseState(2) + measurementState(1) );

  Hinv << cos(poseState(2)+measurementState(1)) , -measurementState(0)*sin(poseState(2)+measurementState(1)) ,
          sin(poseState(2)+measurementState(1)) , measurementState(0)*cos(poseState(2)+measurementState(1));

  covariance = Hinv * measurementUncertainty *Hinv.transpose();
  landmark.set( mean, covariance );

}

double RangeBearingModel::probabilityOfDetection( Pose2d &pose,
						  Landmark2d &landmark,
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
    if( range <= (config.rangeLimMax_ + config.rangeLimBuffer_ ) || 
	range >= (config.rangeLimMin_ - config.rangeLimBuffer_ ) )
      isCloseToSensingLimit = true;
  } 

  return Pd;
}

double RangeBearingModel::clutterIntensity( Measurement2d &z,
					    int nZ){
  return config.uniformClutterIntensity_;
}


double RangeBearingModel::clutterIntensityIntegral( int nZ ){
  double sensingArea_ = 2 * PI * (config.rangeLimMax_ - config.rangeLimMin_);
  return ( config.uniformClutterIntensity_ * sensingArea_ );
}

/************* Implementation of 1d measurement model **************************/

MeasurementModel1d::MeasurementModel1d(){

  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel1d::MeasurementModel1d(Eigen::Matrix<double, 1, 1> &Sr){
  setNoise(Sr);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel1d::MeasurementModel1d(double Sr){
  Measurement1d::Mat S;
  S << Sr;
  setNoise(S);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel1d::~MeasurementModel1d(){}

bool MeasurementModel1d::measure(Pose1d &pose, Landmark1d &landmark, 
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

void MeasurementModel1d::inverseMeasure(Pose1d &pose, Measurement1d &measurement, 
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

double MeasurementModel1d::probabilityOfDetection( Pose1d &pose,
						   Landmark1d &landmark,
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

double MeasurementModel1d::clutterIntensity( Measurement1d &z,
					     int nZ){
  return config.uniformClutterIntensity_;
}


double MeasurementModel1d::clutterIntensityIntegral( int nZ ){
  double sensingLength = config.rangeLimMax_ - config.rangeLimMin_;
  return ( config.uniformClutterIntensity_ * sensingLength );
}
