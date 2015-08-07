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
#include "COLA.hpp"
#include "Particle.hpp"
#include "Pose.hpp"
#include <stdio.h>
#include <string>
#include <vector>

typedef unsigned int uint;

using namespace rfs;

class MM : public Landmark2d{

public:
  
  MM(double x, double y, rfs::TimeStamp t = TimeStamp() ){
    x_(0) = x;
    x_(1) = y;
    this->setTime(t);
  }
  ~MM(){}

  double operator-(const MM& other){

    rfs::Landmark2d::Vec dx = x_ - other.x_;
    return dx.norm();
    
    //return sqrt(mahalanobisDist2( other )); For Mahalanobis distance

  }
  
};


struct COLA_Error{
  
  double error; // combined average error
  double loc; // localization error
  double card; // cardinality error
  TimeStamp t;
  
};


/**
 * \class LogFileReader2dSim
 * \brief A class for reading 2d sim log files and for calculating errors
 */
class LogFileReader2dSim
{
public:

  typedef Particle<Pose2d, GaussianMixture<Landmark2d> > TParticle;
  typedef std::vector< TParticle > TParticleSet;

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
  void readLandmarkGroundtruth(){
    double x, y, t;
    while( fscanf(pGTLandmarkFile, "%lf %lf %lf", &x, &y, &t) == 3){
      map_e_M_.push_back( MM( x, y, TimeStamp(t) ) );
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
	TParticle particle(p_id_, p, p_w_);
	particles_.push_back( particle ); 
	particles_[p_id_].setData( TParticle::PtrData( new GaussianMixture<Landmark2d> ));

	if( fscanf(pParticlePoseFile, "%lf %d %lf %lf %lf %lf", &p_t_, &p_id_, &p_x_, &p_y_, &p_z_, &p_w_) != 6)
	  break;
      }

      // Landmark estimate from highest weighted particle
      double const W_THRESHOLD = 0.75;
      emap_e_M_.clear();
      cardEst_ = 0;
      while(fabs(lmk_t_ - t_currentStep_) < 1e-12){
	if( lmk_pid_ == i_hi_ ){
	  if( lmk_w_ >= W_THRESHOLD ){
	    MM m_e_M(lmk_x_, lmk_y_, lmk_t_);
	    rfs::Landmark2d::Cov mCov;
	    mCov << lmk_sxx_, lmk_sxy_, lmk_sxy_, lmk_syy_;
	    m_e_M.setCov(mCov);
	    emap_e_M_.push_back(m_e_M);
	  }
	  cardEst_ += lmk_w_; 
	}
	
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

    std::vector<MM> map_e_k_M; // observed groundtruth map storage (for Mahananobis distance calculations)
    map_e_k_M.clear();
    for(uint i = 0; i < map_e_M_.size(); i++){
      if( map_e_M_[i].getTime().getTimeAsDouble() <= t_currentStep_ ){
	map_e_k_M.push_back(map_e_M_[i]);
      }
    }
    nLandmarksObservable = map_e_k_M.size();
    
    return cardEst_;
     
   }

  /** Caclculate the error for landmark estimates 
   *  \return COLA errors
   */
  COLA_Error calcLandmarkError(){

    
    double const cutoff = 0.20; // cola cutoff
    double const order = 1.0; // cola order

    std::vector<MM> map_e_k_M; // observed groundtruth map storage (for Mahananobis distance calculations)
    map_e_k_M.clear();
    for(uint i = 0; i < map_e_M_.size(); i++){
      if( map_e_M_[i].getTime().getTimeAsDouble() <= t_currentStep_ ){
	map_e_k_M.push_back(map_e_M_[i]);
      }
    }

    COLA_Error e_cola;
    e_cola.t = t_currentStep_;
    COLA<MM> cola(emap_e_M_, map_e_k_M, cutoff, order);
    e_cola.error = cola.calcError(&(e_cola.loc), &(e_cola.card));

    return e_cola;
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

      poseEst = particles_[i];

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

  TParticleSet particles_;
  std::vector< Pose2d > pose_gt_;
  std::vector<MM> map_e_M_; // groundtruth map storage (for Mahananobis distance calculations)
  std::vector<MM> emap_e_M_; // estimated map storage (for Mahananobis distance calculations)
  double cardEst_; // cardinality estimate

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

    //printf("Time: %f\n", k);
    
    double ex, ey, er, ed;
    
    reader.calcDRError( ex, ey, er, ed);
    fprintf(pDREstErrorFile, "%f   %f   %f   %f   %f\n", k, ex, ey, er, ed);

    reader.calcPoseError( ex, ey, er, ed, false );
    fprintf(pPoseEstErrorFile, "%f   %f   %f   %f   %f\n", k, ex, ey, er, ed);

    //printf("   error x: %f   error y: %f   error rot: %f   error dist: %f\n", ex, ey, er, ed);

    int nLandmarksObserved;
    double cardEst = reader.getCardinalityEst( nLandmarksObserved );
    COLA_Error colaError = reader.calcLandmarkError();
    fprintf(pMapEstErrorFile, "%f   %d   %f   %f\n", k, nLandmarksObserved, cardEst, colaError.error);
 
    //printf("   nLandmarks: %d   nLandmarks estimated: %f   COLA error: %f\n", nLandmarksObserved, cardEst, colaError.loc + colaError.card);
    //printf("--------------------\n");
    
    k = reader.readNextStepData();
  }
  fclose(pPoseEstErrorFile);
  fclose(pMapEstErrorFile);
  fclose(pDREstErrorFile);
  return 0;

}
