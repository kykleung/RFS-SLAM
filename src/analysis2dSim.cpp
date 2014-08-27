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

#include <assert.h> 
#include <boost/filesystem.hpp>
#include "GaussianMixture.hpp"
#include "Landmark.hpp"
#include "HungarianMethod.hpp"
#include "Particle.hpp"
#include "Pose.hpp"
#include <stdio.h>
#include <string>
#include <vector>

using namespace rfs;

/**
 * \class LogFileReader2dSim
 * \brief A class for reading 2d sim log files and for calculating errors
 */
class LogFileReader2dSim
{
public:

  /** Constructor */
  LogFileReader2dSim(const char* logDir){

    std::string filename_gtpose( logDir );
    std::string filename_gtlmk( logDir );
    std::string filename_pose( logDir );
    std::string filename_lmk( logDir );
    std::string filename_dr( logDir );
    
    filename_gtpose += "gtPose.dat";
    filename_gtlmk += "gtLandmark.dat";
    filename_pose += "particlePose.dat";
    filename_lmk += "landmarkEst.dat";
    filename_dr += "deadReckoning.dat";

    pGTPoseFile = fopen(filename_gtpose.data(), "r");
    pGTLandmarkFile = fopen(filename_gtlmk.data(), "r");
    pParticlePoseFile = fopen(filename_pose.data(), "r");
    pLandmarkEstFile = fopen(filename_lmk.data(), "r");
    pDRFile = fopen(filename_dr.data(), "r");

    readLandmarkGroundtruth();

    int nread = fscanf(pParticlePoseFile, "%lf %d %lf %lf %lf %lf", &p_t_, &p_id_, &p_x_, &p_y_, &p_z_, &p_w_);
    //printf("n = %d\n", nread);
    //printf("t = %f   [%d] %f %f %f\n", p_t_, p_id_, p_x_, p_y_, p_w_);

    nread = fscanf(pLandmarkEstFile, "%lf %d %lf %lf %lf %lf %lf %lf\n", 
			   &lmk_t_, &lmk_pid_, &lmk_x_, &lmk_y_, &lmk_sxx_, &lmk_sxy_, &lmk_syy_, &lmk_w_);
    //printf("n = %d\n", nread);
    //printf("t = %f   [%d] %f %f %f\n", lmk_t_, lmk_pid_, lmk_x_, lmk_y_, lmk_w_);
    
  }

  /** Destructor */
  ~LogFileReader2dSim(){

    fclose(pGTPoseFile);
    fclose(pGTLandmarkFile);
    fclose(pParticlePoseFile);
    fclose(pLandmarkEstFile);
    fclose(pDRFile);
    
  }

  /** Read landmark groundtruth data 
   *  \return number of landmarks 
   */
  int readLandmarkGroundtruth(){
    double x, y, t;
    Landmark2d::Vec vm;
    while( fscanf(pGTLandmarkFile, "%lf %lf %lf", &x, &y, &t) == 3){
      vm << x, y;
      map_.push_back( vm );
      mapObsTimestep_.push_back(t);
    }
  }

  /** Read data for the next timestep 
   * \return time for which data was read
   */
  double readNextStepData(){

    if( fscanf(pGTPoseFile, "%lf %lf %lf %lf\n", &t_currentStep_, &rx_, &ry_, &rz_ ) == 4){

      // Dead reckoning
      int rval = fscanf(pDRFile, "%lf %lf %lf %lf\n", &dr_t_, &dr_x_, &dr_y_, &dr_z_);

      // Particles
      particles_.clear();
      particles_.reserve(200);
      w_sum_ = 0;
      double w_hi = 0;
      while(fabs(p_t_ - t_currentStep_) < 1e-12 ){
	// printf("t = %f   [%d] %f %f %f\n", p_t_, p_id_, p_x_, p_y_, p_w_);

	w_sum_ += p_w_;
	if( p_w_ > w_hi ){
	  i_hi_ = p_id_;
	  w_hi = p_w_;
	}

	Pose2d p;
	p[0] = p_x_;
	p[1] = p_y_;
	p[2] = p_z_;
	Particle<Pose2d, GaussianMixture<Landmark2d> > particle(p_id_, p, p_w_);
	particles_.push_back( particle ); 
	particles_[p_id_].setData(new GaussianMixture<Landmark2d>);

	if( fscanf(pParticlePoseFile, "%lf %d %lf %lf %lf %lf", &p_t_, &p_id_, &p_x_, &p_y_, &p_z_, &p_w_) != 6)
	  break;
      }

      // Landmarks
      Landmark2d::Vec vm;
      Landmark2d::Mat Sm;
      while(fabs(lmk_t_ - t_currentStep_) < 1e-12){
	//printf("t = %f   [%d] %f %f %f\n", lmk_t_, lmk_pid_, lmk_x_, lmk_y_, lmk_w_);

	vm << lmk_x_, lmk_y_;
	Sm << lmk_sxx_, lmk_sxy_, lmk_sxy_, lmk_syy_;
	particles_[lmk_pid_].getData()->addGaussian(new Landmark2d(vm, Sm), lmk_w_);

        if(fscanf(pLandmarkEstFile, "%lf %d %lf %lf %lf %lf %lf %lf\n", 
		  &lmk_t_, &lmk_pid_, &lmk_x_, &lmk_y_, &lmk_sxx_, &lmk_sxy_, &lmk_syy_, &lmk_w_) != 8)
	  break;
      }

      return t_currentStep_;
    }
    return -1;
  }

  /** Calculate the cardinality error for landmark estimates
   *  \param[out] nLandmarksObservable the actual number of observed landnmarks up to the current time
   *  \return cardinality estimate
   */
  double getCardinalityEst( int &nLandmarksObservable ){

    std::vector<int> mapObservable;
     for(int n = 0; n < mapObsTimestep_.size(); n++){
       if( mapObsTimestep_[n] <= t_currentStep_ ){
	 mapObservable.push_back(n);
       }
     }
     nLandmarksObservable = mapObservable.size();

     double cardEst = 0;
     for(int i = 0; i < particles_.size(); i++){

       int nGaussians = particles_[i].getData()->getGaussianCount();
       double nLandmarksEst = 0;
       for(int n = 0; n < nGaussians; n++ ){
	 nLandmarksEst += particles_[i].getData()->getWeight(n); 
       }
       cardEst += nLandmarksEst * particles_[i].getWeight() / w_sum_;
     }

     return cardEst;
     
   }

  /** Caclculate the error for landmark estimates 
   *  \return error
   */
  double calcLandmarkError( bool averageError = true){

    double const c2 = 9.0; // cutoff
    double spatialError_i = 0;
    double cardinalityError_i = 0;
    double ospaError = 0;

    std::vector<int> mapObservable;
    for(int n = 0; n < mapObsTimestep_.size(); n++){
      if( mapObsTimestep_[n] <= t_currentStep_ ){
	mapObservable.push_back(n);
      }
    }

    HungarianMethod hm;

    for(int i = 0; i < particles_.size(); i++){

      if( !averageError && i != i_hi_){
	continue;
      }

      double w_i = 1;
      if( averageError ){
	w_i = particles_[i].getWeight();
      }
      
      // For tracking renmaining groundtruth map weight
      int nLandmarks = mapObservable.size();
      double w_remain_gt[ nLandmarks ];
      for(int n = 0; n < nLandmarks; n++ ){
	w_remain_gt[n] = 1;
      }

      // For tracking remaining map estimate weight
      int nGaussians = particles_[i].getData()->getGaussianCount();
      double w_remain_est[ nGaussians ];
      for(int n = 0; n < nGaussians; n++ ){
	w_remain_est[n] = particles_[i].getData()->getWeight(n); 
      }

      // Create matrix of mahalanobis distance between groundtruth and estimated landmark positions
      int E_size = mapObservable.size();
      if(particles_[i].getData()->getGaussianCount() > E_size){
	E_size = particles_[i].getData()->getGaussianCount();
      }
      double** E = new double* [E_size];
      for(int e = 0; e < E_size; e++){
	E[e] = new double[E_size];
	for(int f = 0; f < E_size; f++){
	  E[e][f] = c2;
	}
      }
      for(int n = 0; n < mapObservable.size(); n++){
	for(int m = 0; m < particles_[i].getData()->getGaussianCount(); m++){
	  E[n][m] = particles_[i].getData()->getGaussian(m)->mahalanobisDist2( map_[mapObservable[n]] );
	}
      }
      
      // Allocate memory for a a copy of E, which we can use as we iterate to produce a smaller distance matrix
      double** Er = new double* [E_size];
      for(int e = 0; e < E_size; e++){
	Er[e] = new double[E_size];
      }
      
      std::vector<int> nonZeroEstIdx;
      std::vector<int> nonZeroMapIdx;

      spatialError_i = 0;

      bool smallerThanCutoffDistanceExists = false;
      do{

	nonZeroEstIdx.clear();
	nonZeroMapIdx.clear();
	smallerThanCutoffDistanceExists = false;
	
	for(int n = 0; n < nGaussians; n++ ){
	  if( w_remain_est[n] > 0){
	    nonZeroEstIdx.push_back(n);
	  }
	}
	for(int n = 0; n < nLandmarks; n++ ){
	  if( w_remain_gt[n] > 0){
	    nonZeroMapIdx.push_back(n);
	  }
	}
	int Er_size = nonZeroMapIdx.size();
	if( nonZeroEstIdx.size() > Er_size ){
	  Er_size = nonZeroEstIdx.size();
	}
	if(nonZeroEstIdx.size() == 0 || nonZeroMapIdx.size() == 0){
	  break;
	}
	for(int e = 0; e < Er_size; e++){
	  for(int f = 0; f < Er_size; f++){
	    if(e < nonZeroMapIdx.size() && f < nonZeroEstIdx.size() ){
	      Er[e][f] = E[nonZeroMapIdx[e]][nonZeroEstIdx[f]];
	      if( Er[e][f] > c2 )
		Er[e][f] = c2;
	    }else{
	      Er[e][f] = c2; 
	    }
	  }
	}

	double cost_tmp;
	int match[Er_size];
	bool success = hm.run(Er, Er_size, match, &cost_tmp, false); // Find best linear assignment to minimize distance
	if(!success){
	  printf("particle %d\n", i);
	  for(int e = 0; e < Er_size; e++){
	    for(int f = 0; f < Er_size; f++){
	      printf("%f   ", Er[e][f]);
	    }
	    printf("\n");
	  }
	  printf("\n");
	  return -1;
	}
	
	double err = 0;
	for(int e = 0; e < nonZeroMapIdx.size(); e++){
	  
	  if(match[e] < nonZeroEstIdx.size() && Er[e][match[e]] < c2){

	    smallerThanCutoffDistanceExists = true;

	    double accountedWeight = fmin( w_remain_est[ nonZeroEstIdx[ match[e] ] ] , w_remain_gt[ nonZeroMapIdx[e] ] );
	    w_remain_est[ nonZeroEstIdx[ match[e] ] ] -= accountedWeight;
	    w_remain_gt[ nonZeroMapIdx[e] ] -= accountedWeight;

	    spatialError_i += Er[e][match[e]] * accountedWeight;

	  }

	}

      }while( smallerThanCutoffDistanceExists );
      
      double unaccountedEstWeight = 0;
      double unaccountedMapWeight = 0;
      for(int e = 0; e < nGaussians; e++){
	unaccountedEstWeight += w_remain_est[e];
      }
      for(int e = 0; e < nLandmarks; e++){
	unaccountedMapWeight += w_remain_gt[e];
      }
      cardinalityError_i = fmax(unaccountedMapWeight, unaccountedEstWeight) * c2;
      //printf("Unaccounted map: %f   Unaccounted est: %f\n", unaccountedMapWeight, unaccountedEstWeight);
      //printf("Spatial: %f   Card: %f\n", spatialError_i, cardinalityError_i);

      double ospaError_i = sqrt( (spatialError_i + cardinalityError_i) / nLandmarks ); 
      //printf("OSPA Error: %f   w: %f\n", ospaError, w_i);
      ospaError += (ospaError_i * w_i);

      for(int e = 0; e < E_size; e++ ){
	delete[] E[e];
	delete[] Er[e];
      }
      delete[] E;
      delete[] Er;

    }  

    if( averageError ){
      return (ospaError / w_sum_);
    }else{
      return ospaError;
    }

  }

  /** Calculate the dead reckoning error */
  void calcDRError(double &err_x, double &err_y, double &err_rot, double &err_dist){
    err_x = dr_x_ - rx_;
    err_y = dr_y_ - ry_;
    err_rot = dr_z_ - rz_;
    if(err_rot > PI)
      err_rot -= 2*PI;
    else if(err_rot < -PI)
      err_rot += 2*PI;
    err_dist = sqrt(err_x * err_x + err_y * err_y);
  }


  /** Calculate the error for vehicle pose estimate */
  void calcPoseError(double &err_x, double &err_y, double &err_rot, double &err_dist, bool getAverageError = true){

    err_x = 0;
    err_y = 0;
    err_rot = 0;
    err_dist = 0;

    Pose2d poseEst;
    double w = 1;
    double ex;
    double ey;
    double er;

    for(int i = 0; i < particles_.size(); i++){

      if( !getAverageError && i != i_hi_){
	continue;
      }

      if( getAverageError ){
	w = particles_[i].getWeight();
      }

      particles_[i].getPose(poseEst);

      ex = poseEst[0] - rx_;
      ey = poseEst[1] - ry_;
      er = poseEst[2] - rz_;
      if(er > PI)
	er -= 2 * PI;
      else if(er < -PI)
	er += 2 * PI;

      err_x += ex * w;
      err_y += ey * w;
      err_rot += er * w;
      err_dist += sqrt(ex * ex + ey * ey) * w;

      if( !getAverageError && i == i_hi_){
	break;
      }

    }

    if( getAverageError ){

      err_x /= w_sum_;
      err_y /= w_sum_;
      err_rot /= w_sum_;
      err_dist /= w_sum_;
    }
    
  }

private:

  FILE* pGTPoseFile;       /**< robot pose groundtruth file pointer */
  FILE* pGTLandmarkFile;   /**< landmark groundtruth file pointer */
  FILE* pParticlePoseFile; /**< pose estimate file pointer */
  FILE* pLandmarkEstFile;  /**< landmark estimate file pointer */ 
  FILE* pMapEstErrorFile;  /**< landmark estimate error file pointer */
  FILE* pDRFile;           /**< Dead-reckoning file */

  double mx_; /**< landmark x pos */
  double my_; /**< landmark y pos */

  double t_currentStep_;

  double rx_;
  double ry_;
  double rz_;

  double dr_x_;
  double dr_y_;
  double dr_z_;
  double dr_t_;

  double i_hi_;  /** highest particle weight index */
  double w_sum_; /** particle weight sum */

  std::vector< Particle<Pose2d, GaussianMixture<Landmark2d> > > particles_;
  std::vector< Pose2d > pose_gt_;
  std::vector< Landmark2d::Vec > map_;
  std::vector<double> mapObsTimestep_;

  double p_t_; 
  int p_id_;
  double p_x_; 
  double p_y_; 
  double p_z_; 
  double p_w_;

  double lmk_t_; 
  int lmk_pid_;
  double lmk_x_;
  double lmk_y_;
  double lmk_sxx_;
  double lmk_sxy_;
  double lmk_syy_;
  double lmk_w_;

};

int main(int argc, char* argv[]){

  if( argc != 2 ){
    printf("Usage: analysis2dSim DATA_DIR/\n");
    return 0;
  }
  const char* logDir = argv[1];
  printf("Log directory: %s\n", logDir);

  boost::filesystem::path dir(logDir);
  if(!exists(dir)){
    printf("Log directory %s does not exist\n", logDir);
    return 0;
  }

  std::string filenameLandmarkEstError( logDir );
  std::string filenamePoseEstError( logDir );
  std::string filenameDREstError( logDir );
  filenameLandmarkEstError += "landmarkEstError.dat";
  filenamePoseEstError += "poseEstError.dat";
  filenameDREstError += "deadReckoningError.dat";
  FILE* pMapEstErrorFile = fopen(filenameLandmarkEstError.data(), "w");
  FILE* pPoseEstErrorFile = fopen(filenamePoseEstError.data(), "w");
  FILE* pDREstErrorFile = fopen(filenameDREstError.data(), "w");

  LogFileReader2dSim reader(logDir);

  double k = reader.readNextStepData();
  while( k != -1){

    printf("Time: %f\n", k);
    
    double ex, ey, er, ed;
    
    reader.calcDRError( ex, ey, er, ed);
    fprintf(pDREstErrorFile, "%f   %f   %f   %f   %f\n", k, ex, ey, er, ed);

    reader.calcPoseError( ex, ey, er, ed, false );
    fprintf(pPoseEstErrorFile, "%f   %f   %f   %f   %f\n", k, ex, ey, er, ed);

    printf("   error x: %f   error y: %f   error rot: %f   error dist: %f\n", ex, ey, er, ed);

    int nLandmarksObserved;
    double cardEst = reader.getCardinalityEst( nLandmarksObserved );
    double ospaError = reader.calcLandmarkError( false );
    fprintf(pMapEstErrorFile, "%f   %d   %f   %f\n", k, nLandmarksObserved, cardEst, ospaError);
 
    printf("   nLandmarks: %d   nLandmarks estimated: %f   OSPA error: %f\n", nLandmarksObserved, cardEst, ospaError);
    printf("--------------------\n");
    
    k = reader.readNextStepData();
  }
  fclose(pPoseEstErrorFile);
  fclose(pMapEstErrorFile);
  fclose(pDREstErrorFile);
  return 0;

}
