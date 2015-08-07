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

#include "MeasurementModel_VictoriaPark.hpp"

namespace rfs{

  MeasurementModel_VictoriaPark::MeasurementModel_VictoriaPark(){}

  MeasurementModel_VictoriaPark::~MeasurementModel_VictoriaPark(){}

  MeasurementModel_VictoriaPark::MeasurementModel_VictoriaPark(Eigen::Matrix3d &covZ, double Slb){

    this->setNoise(covZ,Slb);
    config.bufferZonePd_ = 0;
  }

  MeasurementModel_VictoriaPark::MeasurementModel_VictoriaPark(Eigen::Matrix2d &covP, double covR, double Slb){

    Eigen::Matrix3d M;
    M.setZero();
    M.block<2,2>(0,0)=covP;
    M(2,2)=covR;
    this->setNoise(M,Slb);
    config.bufferZonePd_ = 0;
  }

  MeasurementModel_VictoriaPark::MeasurementModel_VictoriaPark(double Sx, double Sy , double SR, double Slb){

    Eigen::Matrix3d M;
    M.setZero();
    M(0,0)=Sx;
    M(1,1)=Sy;
    M(2,2)=SR;
    this->setNoise(M,Slb);  
    config.bufferZonePd_ = 0;
  }

  void MeasurementModel_VictoriaPark::setNoise( Measurement3d::Mat &R , double Slb){
    
    R_=R;
    Slb_=Slb;
    Eigen::Matrix2d cov=R.block<2,2>(0,0);
    rangeBearingModel.setNoise(cov);
  }


  void MeasurementModel_VictoriaPark::inverseMeasure(const Pose2d &pose, const Measurement3d &measurement, Landmark3d &landmark){

    Eigen::Matrix3d covLand,covMeas;
    Eigen::Vector3d meanLand,meanMeas,poseMean;
    TimeStamp t;

    pose.get(poseMean);
    poseMean[2]=poseMean[2]-PI/2;
    Pose2d transformedPose(poseMean, pose.getTime());

    measurement.get(meanMeas,covMeas,t);

    Eigen::Matrix2d cov2d=covMeas.block<2,2>(0,0);
    Eigen::Vector2d mean2d=meanMeas.block<2,1>(0,0);

    Measurement2d meas(mean2d,cov2d,t);
    Landmark2d lan;
    this->rangeBearingModel.inverseMeasure(transformedPose,meas,lan);
    lan.get(mean2d,cov2d,t);
    covLand.setZero();
    covLand.block<2,2>(0,0)=cov2d;
    covLand(2,2)=this->R_(2,2);
    meanLand.setZero();
    meanLand.block<2,1>(0,0)=mean2d;
    meanLand(2)=meanMeas(2);

    landmark.set(meanLand,covLand,t);
  }

  bool MeasurementModel_VictoriaPark::measure( const Pose2d &pose, const Landmark3d &landmark,
					       Measurement3d &measurement, Eigen::Matrix3d *jacobian_wrt_lmk, 
					       Eigen::Matrix3d *jacobian_wrt_pose){

    Eigen::Matrix3d cov,covMeas;
    Eigen::Vector3d mean,meanMeas,poseMean;
    TimeStamp t;

    pose.get(poseMean);
    poseMean[2]=poseMean[2]-PI/2;
    Pose2d transformedPose(poseMean , pose.getTime());

    landmark.get(mean,cov,t);

    Eigen::Matrix2d cov2d=cov.block<2,2>(0,0);
    Eigen::Vector2d mean2d=mean.block<2,1>(0,0);

    Landmark2d lan(mean2d,cov2d);

    Measurement2d meas;
    if(jacobian_wrt_lmk==NULL){
      this->rangeBearingModel.measure(transformedPose,lan,meas,NULL);
      meas.get(mean2d,cov2d,t);
      meanMeas.block<2,1>(0,0)=mean2d;
      meanMeas(2)=mean(2);
      covMeas.setZero();
      covMeas.block<2,2>(0,0)=cov2d;
      covMeas(2,2)=cov(2,2)+this->R_(2,2)+pow(mean2d(0),2)*this->Slb_;

      measurement.set(meanMeas,covMeas,t);

    }else{
      Eigen::Matrix2d jac;
      this->rangeBearingModel.measure(transformedPose,lan,meas,&jac);
      meas.get(mean2d,cov2d,t);
      meanMeas.block<2,1>(0,0)=mean2d;
      meanMeas(2)=mean(2);
      covMeas.setZero();
      covMeas.block<2,2>(0,0)=cov2d;
      covMeas(2,2)=cov(2,2)+this->R_(2,2)+pow(mean2d(0),2)*this->Slb_;
      jacobian_wrt_lmk->setZero();
      jacobian_wrt_lmk->block<2,2>(0,0)=jac;
      (*jacobian_wrt_lmk)(2,2)=1;

      measurement.set(meanMeas,covMeas,t);
    }
    return true;
  }

  double MeasurementModel_VictoriaPark::probabilityOfDetection( const Pose2d &pose,
								const Landmark3d &landmark,
								bool &isCloseToSensingLimit){


    Eigen::Vector3d circleMean,measMean,circleMean2,poseMean;
    Eigen::Matrix3d circleCov;
    landmark.get(circleMean,circleCov);
    pose.get(poseMean);
    Measurement3d measurement;
    this->measure(pose,landmark,measurement);
    measurement.get(measMean);
    double distancetocircle = sqrt(pow(measMean[0],2)+pow(measMean[1],2));
    double angle = atan2(measMean[1],measMean[0])+poseMean[2];

    Eigen::Vector2d perp_direction ,direction;
    perp_direction << -sin(angle) , cos(angle);
    direction << cos(angle) , sin(angle);
    double std = perp_direction.transpose()*circleCov.block<2,2>(0,0)*perp_direction;
    std= 3*sqrt(std);
    std=std::max(std,0.2);
    circleMean2(2)=circleMean(2);
    circleMean2.block<2,1>(0,0) = circleMean.block<2,1>(0,0)+ std*perp_direction;

    Landmark3d temp_landmark(circleMean2,circleCov);

    double Pd0,Pd1,Pd2;
    std::vector<double> pd_vect;

    for (int i=1;(i-1)*(2*circleMean(2))<std;i++){
      circleMean2.block<2,1>(0,0) = circleMean.block<2,1>(0,0)+i*2*circleMean(2)*perp_direction;
      temp_landmark.set(circleMean2,circleCov);
      pd_vect.push_back(probabilityOfDetection2(pose,temp_landmark,isCloseToSensingLimit));
      circleMean2.block<2,1>(0,0) = circleMean.block<2,1>(0,0)-i*2*circleMean(2)*perp_direction;
      temp_landmark.set(circleMean2,circleCov);
      pd_vect.push_back(probabilityOfDetection2(pose,temp_landmark,isCloseToSensingLimit));	  	
    }
    pd_vect.push_back(probabilityOfDetection2(pose,landmark,isCloseToSensingLimit));


    if(*std::min_element(pd_vect.begin(),pd_vect.end())==0 && *std::max_element(pd_vect.begin(),pd_vect.end())>0){

      isCloseToSensingLimit= true;
    }

    return *std::max_element(pd_vect.begin(),pd_vect.end());
  }


  double MeasurementModel_VictoriaPark::probabilityOfDetection2( const Pose2d &pose,
								 const Landmark3d &landmark,
								 bool &isCloseToSensingLimit){
    //std::cout<< "pD"<<std::endl;
    isCloseToSensingLimit=false;
    Eigen::Vector3d circleMean,measMean;
    Eigen::Matrix3d circleCov;
    landmark.get(circleMean,circleCov);

    Measurement3d measurement;

    this->measure(pose,landmark,measurement);
    measurement.get(measMean);
    double distancetocircle = measMean[0];
    double angle = measMean[1];
  
    if (angle > config.bearingLimitMax_ || angle < config.bearingLimitMin_ || 
	distancetocircle < config.rangeLimMin_ || distancetocircle > config.rangeLimMax_){
      return 0;
    }
    // return config.probabilityOfDetection_[config.probabilityOfDetection_.size()-1];

    double modified_radius=measMean[2]/2;//+sqrt(circleCov(0,0)+circleCov(1,1)+circleCov(2,2));
    double gamma=atan(modified_radius/measMean[0]);

    int maxNumPoints =(int)floor(2*gamma*720.0/(2*PI));

    if(config.probabilityOfDetection_.size()>maxNumPoints && config.probabilityOfDetection_[maxNumPoints]==0){
      return 0;
    }
    if(config.probabilityOfDetection_.size()>maxNumPoints && config.probabilityOfDetection_[maxNumPoints] < config.bufferZonePd_){
      isCloseToSensingLimit =true;
    }

    int minb=(int)ceil((angle-gamma)*720.0/(2*PI));
    int maxb=minb+maxNumPoints;

    while(minb>=720)minb-=720;
    while(minb<0)minb+=720;
    while(maxb>=720)maxb-=720;
    while(maxb<0)maxb+=720;

    int numPoints=0;

    double minrange=distancetocircle - modified_radius - 6*0.03;// - sqrt(std::max(circleCov(0,0),circleCov(1,1))+circleCov(2,2));

    if((maxb-minb+720)%720>0){
      for(int bearing=minb;bearing!=maxb;bearing=(bearing+1)%720){

	if(this->laserscan_[bearing]>minrange || this->laserscan_[bearing]==0)
	  numPoints++;
      }
    }

    int numPoints_small=numPoints;

    if(numPoints>=config.probabilityOfDetection_.size()){
      numPoints=config.probabilityOfDetection_.size()-1;
    }

    if(config.probabilityOfDetection_[numPoints]==0)
      isCloseToSensingLimit = false;
    return config.probabilityOfDetection_[numPoints];
  }

  void MeasurementModel_VictoriaPark::setLaserScan(const std::vector<double> &laserscan){

    laserscan_=laserscan;

    double FoVArea=0;
    for(int i=1;i<laserscan_.size();i++){
     
      FoVArea+=laserscan_[i]*laserscan_[i-1];
    }
    FoVArea+=laserscan_[0]*laserscan_[laserscan_.size()-1];
    FoVArea*=sin(PI/360)/2;
     
    clutterIntensity_=config.expectedClutterNumber_/FoVArea;
     
  }

  double MeasurementModel_VictoriaPark::clutterIntensity( Measurement3d &z, int nZ ){

    return clutterIntensity_;
  }

  double MeasurementModel_VictoriaPark::clutterIntensityIntegral( int nZ){

    return config.expectedClutterNumber_;
  }

}
