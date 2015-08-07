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

#include "MeasurementModel_XY.hpp"

namespace rfs
{

MeasurementModel_XY::MeasurementModel_XY(){

  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}


MeasurementModel_XY::MeasurementModel_XY(::Eigen::Matrix2d &covZ){

  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_XY::MeasurementModel_XY(double Sx, double Sy){

  Eigen::Matrix2d covZ;
  covZ <<  Sx, 0, 0, Sy;
  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_XY::~MeasurementModel_XY(){}

bool MeasurementModel_XY::measure(const Pose2d &pose, 
				      const Landmark2d &landmark, 
				      Measurement2d &measurement,
				      Eigen::Matrix2d *jacobian_wrt_lmk,
				      Eigen::Matrix<double, 2, 3> *jacobian_wrt_pose){
  
  Eigen::Vector3d robotPose;
  Eigen::Vector2d mean, landmarkState;
  Eigen::Matrix2d H_lmk, landmarkUncertainty, cov;
  Eigen::Matrix<double, 2, 3> H_robot;
  Eigen::Matrix3d robotUncertainty;
  double dx_I, dy_I, dx_R, dy_R, c_RI, s_RI, range, range2, bearing;

  pose.get(robotPose, robotUncertainty);
  landmark.get(landmarkState,landmarkUncertainty);

  dx_I = landmarkState(0) - robotPose(0);
  dy_I = landmarkState(1) - robotPose(1);
  c_RI = cos(robotPose(2));
  s_RI = sin(robotPose(2));

  range2 = pow(dx_I, 2) + pow(dy_I, 2) ;
  range = sqrt(range2);

  dx_R =  c_RI * dx_I + s_RI * dy_I;
  dy_R = -s_RI * dx_I + c_RI * dy_I;
  mean << dx_R, dy_R ;

  H_lmk <<  c_RI, s_RI,
           -s_RI, c_RI;

  H_robot << -c_RI, -s_RI, -dx_I*s_RI + dy_I*c_RI,
              s_RI, -c_RI, -dx_I*c_RI - dy_I*s_RI;

  cov = H_lmk * landmarkUncertainty * H_lmk.transpose() + H_robot * robotUncertainty * H_robot.transpose() + R_;
  measurement.set(mean, cov);

  if(jacobian_wrt_lmk != NULL)
    *jacobian_wrt_lmk = H_lmk;

  if(jacobian_wrt_pose != NULL)
    *jacobian_wrt_pose = H_robot;

  if(range > config.rangeLimMax_ || range < config.rangeLimMin_)
    return false;
  else
    return true;
}

void MeasurementModel_XY::inverseMeasure(const Pose2d &pose, 
					 const Measurement2d &measurement, 
					 Landmark2d &landmark){
  Eigen::Vector3d poseState;
  Eigen::Vector2d measurementState, mean;
  Eigen::Matrix2d measurementUncertainty, covariance, Hinv;
 
  pose.get(poseState);
  measurement.get(measurementState);
  this->getNoise(measurementUncertainty); 

  double c_RI = cos(poseState(2));
  double s_RI = sin(poseState(2));

  mean << poseState(0) + c_RI * measurementState(0) - s_RI * measurementState(1),
          poseState(1) + s_RI * measurementState(0) + c_RI * measurementState(1);

  Hinv << c_RI, -s_RI,
          s_RI, c_RI;
  
  covariance = Hinv * measurementUncertainty * Hinv.transpose();
  landmark.set( mean, covariance );

}

double MeasurementModel_XY::probabilityOfDetection( const Pose2d &pose,
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

double MeasurementModel_XY::clutterIntensity( Measurement2d &z,
					    int nZ){
  return config.uniformClutterIntensity_;
}


double MeasurementModel_XY::clutterIntensityIntegral( int nZ ){
  double sensingArea_ = 2 * PI * (config.rangeLimMax_ - config.rangeLimMin_);
  return ( config.uniformClutterIntensity_ * sensingArea_ );
}

}
