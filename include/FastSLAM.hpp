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

#ifndef FASTSLAM_HPP
#define FASTSLAM_HPP

#include "Timer.hpp"
#include <Eigen/Core>
#include "GaussianMixture.hpp"
#include "HungarianMethod.hpp"
#include "KalmanFilter.hpp"
#include "MurtyAlgorithm.hpp"
#include "ParticleFilter.hpp"
#include <math.h>
#include <vector>
#include <stdio.h>
#include <list>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rfs{

/**
 *  \class FastSLAM
 *  \brief Factored Solution to SLAM
 *  
 *  This is an implementation of the FastSLAM v1.0 algorithm and the MH-FastSLAM algorithm
 *
 *  INPROCEEDINGS{Montemerlo02a,
 *  AUTHOR         = {Montemerlo, M. and Thrun, S. and Koller, D. and 
 *                   Wegbreit, B.},
 *  TITLE          = {{FastSLAM}: {A} Factored Solution to the Simultaneous 
 *                   Localization and Mapping Problem},
 *  YEAR           = {2002},
 *  BOOKTITLE      = {Proceedings of the AAAI National Conference on 
 *                   Artificial Intelligence},
 *  PUBLISHER      = {AAAI},
 *  ADDRESS        = {Edmonton, Canada}
 *  }
 *
 *  \tparam RobotProcessModel A robot process model derived from ProcessModel
 *  \tparam LmkProcessModel A landmark process model derived from ProcessModel
 *  \tparam MeasurementModel A sensor model derived from MeasurementModel
 *  \tparam KalmanFilter A Kalman filter that uses LmkProcessModel and MeasurementModel
 *  \author Keith Leung, Felipe Inostroza
 */

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
class FastSLAM : public ParticleFilter<RobotProcessModel, MeasurementModel, 
				       GaussianMixture< typename MeasurementModel::TLandmark > >
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef typename RobotProcessModel::TState TPose;
  typedef Trajectory<TPose> TTrajectory;
  typedef typename RobotProcessModel::TInput TInput;
  typedef typename MeasurementModel::TLandmark TLandmark;
  typedef typename MeasurementModel::TMeasurement TMeasurement;
  typedef GaussianMixture<TLandmark> TGM;
  typedef typename TGM::Gaussian TGaussian;
  typedef ParticleFilter<RobotProcessModel, MeasurementModel, TGM > TPF;
  typedef typename TPF::TParticle TParticle;
  typedef typename TPF::TParticleSet TParticleSet;

  /**
   * \class LandmarkCandidate
   * \brief New landmark candidate
   */
  class LandmarkCandidate : public TLandmark{
  public:
    uint nSupportingMeasurements;
    uint nChecks;
  };

  /** 
   * \brief Configurations for this RBPHDFilter 
   */
  struct Config{

    /** Minimum updates betwen resampling of particles*/
    int minUpdatesBeforeResample_;

    /** Minimum numnber of measurements before resampling of particles*/
    int minMeasurementsBeforeResample_;

    /** If true, timing information is written to the console every update*/
    bool reportTimingInfo_;

    /** The probability of existence that is initially assigned to a new landmark */
    double landmarkExistencePrior_;

    /** The log odds threshold for eliminating a landmark from the map*/
    double mapExistencePruneThreshold_;

    /** Minimum log measurement likelihood for numerical stability */
    double minLogMeasurementLikelihood_;

    /** Maximum number of particles */
    int nParticlesMax_;

    /** Maximum number of data associations per particle for MH-FastSLAM**/
    unsigned int maxNDataAssocHypotheses_;

    /** Maximum measurement likelihood difference to the previous best data association hypothesis for a given association to be allowed to spawn a new particle in MH-FastSLAM */
    double maxDataAssocLogLikelihoodDiff_;

    /**  Mahalanobis distance threshold less than which a measurement counts towards supporting a new landmark candidate */
    double landmarkCandidateMeasurementSupportDist_;

    /**  Number of supporting measurements required to generate a new landmark */
    uint landmarkCandidateMeasurementCountThreshold_;

    /**  Number of latest observations below or equal to which birth Gaussians are added regardless or measurement evidence.
     *   This is to accomodate for the case where the robot enters an area sparesly populated with landmarks.
     */
    uint landmarkCandidateCurrentMeasurementCountThreshold_;

    /**  Number of checks before candidate birth Gaussian is removed if there are not enough supporting measurements */
    uint landmarkCandidateMeasurementCheckThreshold_;

    /** Weight above which a landmark is considered permanent */
    double landmarkLockWeight_;

    /** Map pruning occurs if number of measurements is equal or above this threshold during update */
    uint pruningMeasurementsThreshold_;

  } config;

  
  /**
   * \brief Elapsed timing information 
   */
  struct TimingInfo {
    long long predict_wall;
    long long predict_cpu;
    long long mapUpdate_wall;
    long long mapUpdate_cpu;
    long long dataAssoc_wall;
    long long dataAssoc_cpu;
    long long mapUpdate_KF_wall;
    long long mapUpdate_KF_cpu;
    long long weighting_wall;
    long long weighting_cpu;
    long long mapManage_wall;
    long long mapManage_cpu;
    long long particleResample_wall;
    long long particleResample_cpu;
  } timingInfo_;

  /** 
   * Constructor 
   * \param n number of particles
   */
  FastSLAM(int n);

  /** Destructor */
  ~FastSLAM();

  /** 
   * Get the landmark process model
   * \return pointer to the landmark process model
   */
  LmkProcessModel* getLmkProcessModel();

  /**
   * Predict the robot trajectory using the lastest odometry data
   * \param[in] u input 
   * \param[in] dT size of timestep;
   * \param[in] useModelNoise use the additive noise for the motion model
   * \param[in] useInputNoise use the noise for the inputs
   */
  void predict( TInput u, const TimeStamp &dT,
		bool useModelNoise = true,
		bool useInputNoise = false);

  /**
   * Update the map, calculate importance weighting, and perform resampling if necessary
   * \param[in] Z set of measurements to use for the update, placed in a std vector, which
   * gets cleared after the function call. 
   */
  void update( std::vector<TMeasurement> &Z);

  /**
   * Get the size of the Gaussian mixture for a particle
   * \param[in] i particle index
   * \return size if index is valid, else -1
   */
  int getGMSize(int i);

  /**
   * Get the position, covariance, and weight of a Gaussian in particle i's Gaussin mixture
   * \param[in] i particle index
   * \param[in] m Gaussian index
   * \param[out] u mean
   * \param[out] S covariance
   * \param[out] w log-odds of existence
   * \return false if the indices specified are invalid 
   */ 
  bool getLandmark(const int i, const int m, 
		   typename TLandmark::Vec &u,
		   typename TLandmark::Mat &S,
		   double &w);
  
  /**
   * Get the pointer to the Kalman Filter used for updating the map
   * \return pointer to the Kalman Filter
   */
  KalmanFilter* getKalmanFilter();

  /** Function for initiating particles during startup 
   *  \param[in] i particle index
   *  \param[in] p particle pose
   */
  void setParticlePose(int i, TPose &p);

  /** 
   * Get elapsed timing information for various steps of the filter 
   */
  TimingInfo* getTimingInfo();


private:

  std::vector< std::list<LandmarkCandidate> > landmarkCandidates_;
  std::vector< uint > nLandmarksInFOV_;

  uint nThreads_;

  std::vector<KalmanFilter> kfs_; /**< Kalman filters (one for each thread)*/
  LmkProcessModel *lmkModelPtr_; /**< pointer to landmark process model */

  int nParticles_init_;

  bool resampleOccured_;

  Timer timer_predict_; /**< Timer for prediction step */
  Timer timer_mapUpdate_; /**< Timer for map update */
  Timer timer_dataAssoc_; /**< Timer for data assoication */
  Timer timer_mapUpdate_KF_; /**< Timer for map update with the Kalman Filter */
  Timer timer_weighting_; /**< Timer for particle weighting */
  Timer timer_mapManage_; /**< Timer for map management */
  Timer timer_particleResample_; /**<Timer for particle resampling */ 

  unsigned int nUpdatesSinceResample_; /**< Number of updates performed since the last resmaple */
  unsigned int nMeasurementsSinceResample_; /**< Number of measurements processed since the last resample */

  /**
   * Update the map with the measurements in measurements_
   * Existing landmarks with probability of detection > 0 will have their Gaussian
   * mixture weight reduced to account for missed detection.
   * For every landmark-measurement pair with probability of detection > 0,
   * a new landmark will be created. 
   * /param[in] particleIdx index of the particle for which the map is to be updated
   */
  void updateMap(const uint particleIdx);

  /**
   * Resample the particles, along with their individual maps,  according to their 
   * importance weightings.
   */
  void resampleWithMapCopy();

  /** 
   * Importance weighting. Overrides the abstract function in ParticleFilter
   * \note For this FastSLAM algorithm, we perform weighting as part of 
   * mapUpdate() to be a little more efficient. Therefore, this function is 
   * not called at all. However, we still need to overwrite the virtual function.
   * \param[in] idx particle index
   */
  void importanceWeighting(const uint idx){}

};

////////// Implementation //////////

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::FastSLAM(int n)
: ParticleFilter<RobotProcessModel, MeasurementModel, 
		 GaussianMixture< typename MeasurementModel::TLandmark > >(n)
{

  nThreads_= 1;

  #ifdef _OPENMP
  nThreads_ = omp_get_max_threads();
  #endif

  lmkModelPtr_ = new LmkProcessModel;
  kfs_ = std::vector<KalmanFilter>(nThreads_, 
				   KalmanFilter(lmkModelPtr_, this->getMeasurementModel() ) );

  nParticles_init_ = n;
  
  for(int i = 0; i < n; i++){
    this->particleSet_[i]->setData( typename TParticle::PtrData( new GaussianMixture<TLandmark>() ) );
  }

  config.minUpdatesBeforeResample_ = 1;
  config.minMeasurementsBeforeResample_ = 1;
  config.reportTimingInfo_ = false;
  config.landmarkExistencePrior_ = 0.5;
  config.mapExistencePruneThreshold_ = -3.0;
  config.minLogMeasurementLikelihood_ = -10.0;
  config.nParticlesMax_ = n * 3;
  config.maxNDataAssocHypotheses_ = 1;
  config.maxDataAssocLogLikelihoodDiff_ = 5;
  config.landmarkCandidateMeasurementSupportDist_ = 1;
  config.landmarkCandidateMeasurementCountThreshold_ = 1;
  config.landmarkCandidateCurrentMeasurementCountThreshold_ = 1;
  config.landmarkCandidateMeasurementCheckThreshold_ = 2;
  config.landmarkLockWeight_ = 10;
  config.pruningMeasurementsThreshold_ = 0;

  nUpdatesSinceResample_ = 0;
  nMeasurementsSinceResample_ = 0;

  landmarkCandidates_.resize(n);
  nLandmarksInFOV_.resize(n);
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::~FastSLAM(){
  for(int i = 0; i < this->nParticles_; i++){
    this->particleSet_[i]->deleteData(); //delete maps_[i];
  }
  delete lmkModelPtr_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
LmkProcessModel* FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getLmkProcessModel(){
  return lmkModelPtr_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::predict( TInput u, 
											      const TimeStamp &dT,
											      bool useModelNoise,
											      bool useInputNoise){

  timer_predict_.resume();

  // propagate particles
  this->propagate(u, dT, useModelNoise, useInputNoise, true); // true for keeping trajecotory 

  // propagate landmarks
  for( int i = 0; i < this->nParticles_; i++ ){
    for( int m = 0; m < this->particleSet_[i]->getData()->getGaussianCount(); m++){
      TLandmark *plm;
      this->particleSet_[i]->getData()->getGaussian(m, plm);
      lmkModelPtr_->staticStep(*plm, *plm, dT);
    }
  }

  timer_predict_.stop();
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::update( std::vector<TMeasurement> &Z){

  const unsigned int startIdx = 0;
  const unsigned int stopIdx = this->nParticles_;

  if(landmarkCandidates_.size() < this->nParticles_ * config.maxNDataAssocHypotheses_)
    landmarkCandidates_.resize(this->nParticles_ * config.maxNDataAssocHypotheses_);
  if(nLandmarksInFOV_.size() < this->nParticles_ * config.maxNDataAssocHypotheses_)
    nLandmarksInFOV_.resize(this->nParticles_ * config.maxNDataAssocHypotheses_);

  nUpdatesSinceResample_++;

  this->setMeasurements( Z ); // Z gets cleared after this call, measurements now stored in this->measurements_
  if(this->measurements_.size() == 0)
    return;
  nMeasurementsSinceResample_ += this->measurements_.size();

  // make sure any setting changes to the Kalman Filter are set for all threads
  if(nThreads_>1){
    for(int j = 1; j < nThreads_; j++){
      kfs_[j]=kfs_[0];
    }
  }

  ////////// Map Update and Particle Weighting//////////
  timer_mapUpdate_.resume();

  #pragma omp parallel
  {
    #pragma omp for
    for(unsigned int i = startIdx; i < stopIdx; i++){
      updateMap(i);
    }
  }

  timer_mapUpdate_.stop();

  //////////// Particle resampling //////////
  timer_particleResample_.resume();
  resampleWithMapCopy();
  timer_particleResample_.stop();

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::updateMap(const uint particleIdx){

  int threadnum = 0;
  #ifdef _OPENMP
    threadnum = omp_get_thread_num();
  #endif

  const unsigned int nZ = this->measurements_.size();
  const uint i = particleIdx;

  //----------  1. Data Association --------------------
  #ifndef _OPENMP
  timer_dataAssoc_.resume();
  #endif
  
  // Look for landmarks within sensor range
  const TPose *pose = this->particleSet_[i]->getPose();
  unsigned int nM = this->particleSet_[i]->getData()->getGaussianCount();
  std::vector<int> idx_inRange;
  std::vector<double> pd_inRange;
  std::vector<TLandmark*> lm_inRange;
  for( int m = 0; m < nM; m++ ){
    TLandmark* lm = this->particleSet_[i]->getData()->getGaussian(m);
    bool closeToLimit = false;
    double pd = this->pMeasurementModel_->probabilityOfDetection(*pose, *lm, closeToLimit); 
    if( pd != 0 || closeToLimit ){
      idx_inRange.push_back(m);
      lm_inRange.push_back(lm);
      pd_inRange.push_back(pd);
    }
  }
  nM = lm_inRange.size();

  // Loglikelihood table for Data Association
  unsigned int nMZ = nM;
  if(nZ > nM){
    nMZ = nZ;
  }
  double** likelihoodTable;
  CostMatrix likelihoodMat(likelihoodTable, nMZ);
  for( int m = 0; m < nMZ; m++ ){
    for(int z = 0; z < nMZ; z++){
      likelihoodTable[m][z] = config.minLogMeasurementLikelihood_;
    }
  }

  // Fill in table for landmarks within range    
  for(unsigned int m = 0; m < nM; m++){
    TLandmark* lm = lm_inRange[m]; // landmark position estimate
    TMeasurement measurement_exp; // expected measurement      
    bool isValidExpectedMeasurement = this->pMeasurementModel_->measure( *pose , *lm , measurement_exp);
    for(int z = 0; z < nZ; z++){
      if( isValidExpectedMeasurement ){
	likelihoodTable[m][z] = fmax(config.minLogMeasurementLikelihood_, 
				     log(measurement_exp.evalGaussianLikelihood(this->measurements_[z])));
      }
    }
  }
 
  likelihoodMat.reduce(config.minLogMeasurementLikelihood_);
  double** likelihoodTableReduced;
  int* assignments_fixed = new int[nMZ];
  double score_reduced;
  int* mRemap;
  int* zRemap;
  int nMZReduced = likelihoodMat.getCostMatrixReduced(likelihoodTableReduced, assignments_fixed, &score_reduced, mRemap, zRemap);
   
  // Use Hungaian Method and Murty's Method for k-best data association
  unsigned int nH = 0; // number of data association hypotheses (we will create a new particle for each)
  double logLikelihoodSum = 0;
  Murty murty(likelihoodTableReduced, nMZReduced);
  std::vector<int*> da; // data association hypotheses
  if( nMZReduced == 0 ){
    int* da_current = new int[nMZ];
    for(unsigned int m = 0; m < nM; m++){
      da_current[m] = assignments_fixed[m];
    }
    da.push_back( da_current );
  }else{
    while(nH < config.maxNDataAssocHypotheses_){

      Murty::Assignment daVar;
      unsigned int nH_old = nH;
      nH = murty.findNextBest(daVar, logLikelihoodSum);
      if(nH == -1){
	nH = nH_old;
	break;
      }
      if(murty.getBestScore() - logLikelihoodSum >= config.maxDataAssocLogLikelihoodDiff_){
	nH--;
	break;
      }

      int* da_current = new int[nMZ];
      for(unsigned int m = 0; m < nM; m++){
	da_current[m] = assignments_fixed[m];
      }
      for(unsigned int m = 0; m < nMZReduced; m++){
	int z_o = zRemap[ daVar[m] ];
	int m_o = mRemap[m];
	if(z_o < nZ){
	  da_current[ m_o ] = z_o;
	}else{
	  da_current[ m_o ] = -2;
	}
      }
      da.push_back(da_current);

    }
  }

  delete[] assignments_fixed;
  
  // particle indices for update
  unsigned int pi[nH];
  pi[0] = i; 
  if(nH > 1){
    double newWeight = this->particleSet_[i]->getWeight() / nH;
    this->particleSet_[i]->setWeight(newWeight); 
    #pragma omp critical(increaseParticles)
    {
      this->copyParticle(i, nH-1, newWeight);
      for(unsigned int h = 1; h < nH; h++){
	pi[h] = this->nParticles_ - h;
	if(resampleOccured_){
	   landmarkCandidates_[pi[h]] = landmarkCandidates_[pi[0]];
	}
      }
    }
  }

  #ifndef _OPENMP
  timer_dataAssoc_.stop();
  #endif
  
  for(unsigned int h = 0; h < nH; h++){
    
    //----------  2. Kalman Filter map update ----------

    double nExpectedClutter = this->pMeasurementModel_->clutterIntensityIntegral(nZ);
    double probFalseAlarm = nExpectedClutter / nZ;
    double p_exist_given_Z = 0;
    double logParticleWeight = 0;

    nLandmarksInFOV_[pi[h]] = 0;
    bool zUsed[nZ];
    for(unsigned int z = 0; z < nZ; z++){
      zUsed[z] = false;
    }

    #ifndef _OPENMP
    timer_mapUpdate_KF_.resume();
    #endif

    for(unsigned int m = 0; m < nM; m++){
      
      TLandmark* lm = this->particleSet_[pi[h]]->getData()->getGaussian( idx_inRange[m] );

      // Update landmark estimate m with the associated measurement
      int z = da[h][m];
      bool isUpdatePerformed = false;
      if(z < nZ && z >= 0 && likelihoodTable[m][z] > config.minLogMeasurementLikelihood_){
	isUpdatePerformed = kfs_[threadnum].correct(*pose, this->measurements_[z], *lm, *lm);
      }
      
      // calculate change to existence probability      
      double w = this->particleSet_[pi[h]]->getData()->getWeight( idx_inRange[m] );
      if(isUpdatePerformed){

	nLandmarksInFOV_[pi[h]]++;
	zUsed[z] = true; // This flag is for new landmark creation
	logParticleWeight += likelihoodTable[m][z];
	p_exist_given_Z = ((1 - pd_inRange[m]) * probFalseAlarm * config.landmarkExistencePrior_ + pd_inRange[m] * config.landmarkExistencePrior_) /
	  (probFalseAlarm + (1 - probFalseAlarm) * pd_inRange[m] * config.landmarkExistencePrior_); 

      }else{ // landmark estimate m not updated
	  
	p_exist_given_Z = ((1 - pd_inRange[m]) * config.landmarkExistencePrior_) /
	  ((1 - config.landmarkExistencePrior_) + (1 - pd_inRange[m]) * config.landmarkExistencePrior_);

	if(w > config.landmarkLockWeight_)
	  p_exist_given_Z = 0.5;

      }
      
      w += log( (p_exist_given_Z) / (1 - p_exist_given_Z) ); 
      this->particleSet_[pi[h]]->getData()->setWeight(idx_inRange[m], w);
    }
     
    #ifndef _OPENMP
    timer_mapUpdate_KF_.stop();
    timer_mapManage_.resume();
    #endif
  
    //---------- 3. Map Management (Add and remove landmarks)  ------------
    //if(nLandmarksInFOV_[pi[h]] >= 3 && nZ >= 3)
    if(nZ >= config.pruningMeasurementsThreshold_)
      int nRemoved = this->particleSet_[pi[h]]->getData()->prune(config.mapExistencePruneThreshold_); 
    
    double newLandmarkWeight = log(config.landmarkExistencePrior_ / (1 - config.landmarkExistencePrior_));

    for(unsigned int z = 0; z < nZ; z++){
      if(!zUsed[z]){ // Create new landmarks with inverse measurement model with unused measurements

	// get measurement
	TMeasurement unused_z = this->measurements_[z];

	// check to see if measurement corresponds closely with a Landmark in the landmark candidate list
	bool isNewCandidate = true;
      
	for( typename std::list<LandmarkCandidate>::iterator it = landmarkCandidates_[pi[h]].begin();
	     it != landmarkCandidates_[pi[h]].end(); it++ ){
	  
	  TMeasurement z_exp;
	  TPose x = *(this->particleSet_[pi[h]]);
	  this->pMeasurementModel_->measure(x, *it, z_exp);
	  double d2 = z_exp.mahalanobisDist2( unused_z );
	  if(d2 <= config.landmarkCandidateMeasurementSupportDist_ * config.landmarkCandidateMeasurementSupportDist_){
	    kfs_[0].correct(x, unused_z, *it, *it);
	    (it->nSupportingMeasurements)++;
	    isNewCandidate = false;
	    break;
	  }
	}

	if(isNewCandidate){

	  // use inverse measurement model to get landmark
	  TPose robot_pose = *(this->particleSet_[pi[h]]);
	  LandmarkCandidate c;
	  c.nSupportingMeasurements = 1;
	  c.nChecks = 0;
	  this->pMeasurementModel_->inverseMeasure( robot_pose, unused_z, c );

	  if(config.landmarkCandidateMeasurementCountThreshold_ == 1 ||
	     nLandmarksInFOV_[pi[h]] <= config.landmarkCandidateCurrentMeasurementCountThreshold_){
	    // add birth landmark to Gaussian mixture (last param = true to allocate mem)
	    this->particleSet_[pi[h]]->getData()->addGaussian( &c, newLandmarkWeight, true);	  
	  }else{
	    landmarkCandidates_[pi[h]].push_back(c);
	  }
	
	}  
	
	// Check through candidate list to see if any candidates should be made into a real birth Gaussian
	for( typename std::list<LandmarkCandidate>::iterator it = landmarkCandidates_[pi[h]].begin();
	     it != landmarkCandidates_[pi[h]].end(); it++ ){
	  it->nChecks++;
	  while( it->nSupportingMeasurements >= config.landmarkCandidateMeasurementCountThreshold_ || 
		 it->nChecks > config.landmarkCandidateMeasurementCheckThreshold_ ||
		 nLandmarksInFOV_[pi[h]] <= config.landmarkCandidateCurrentMeasurementCountThreshold_){
	    
	    if (it->nSupportingMeasurements >= config.landmarkCandidateMeasurementCountThreshold_){
	      
	      this->particleSet_[pi[h]]->getData()->addGaussian( &(*it), newLandmarkWeight * it->nChecks, true);

	    }else if( nLandmarksInFOV_[pi[h]] <= config.landmarkCandidateCurrentMeasurementCountThreshold_){
	      this->particleSet_[pi[h]]->getData()->addGaussian( &(*it), newLandmarkWeight * it->nChecks, true);
	    }
	    it = landmarkCandidates_[pi[h]].erase( it );
	    if( it != landmarkCandidates_[pi[h]].end() ){
	      it->nChecks++;
	    }else{
	      break;
	    }
	  }
	}
	
      }
    }

    #ifndef _OPENMP
    timer_mapManage_.stop();
    timer_weighting_.resume();
    #endif

    //---------- 4. Importance Weighting --------------
    // Some of the work has already been done in the map update step
    double prev_p_weight = this->particleSet_[pi[h]]->getWeight();
    this->particleSet_[pi[h]]->setWeight( prev_p_weight * exp(logParticleWeight) );

    #ifndef _OPENMP
    timer_weighting_.stop();
    #endif

  }

  //---------- 5. Cleanup - Free memory ---------
  for( int d = 0; d < da.size(); d++ ){
    delete[] da[d];
  }

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::resampleWithMapCopy(){

  resampleOccured_ = false;

  if( this->nParticles_ > config.nParticlesMax_){
    resampleOccured_ = this->resample( nParticles_init_, true );
  }else if( nUpdatesSinceResample_ >= config.minUpdatesBeforeResample_ && 
	    nMeasurementsSinceResample_ >= config.minMeasurementsBeforeResample_){
    resampleOccured_ = this->resample( nParticles_init_ );
  }

  if( resampleOccured_ ){
    nUpdatesSinceResample_ = 0;
    nMeasurementsSinceResample_ = 0;
  }else{
    this->normalizeWeights();
  }

  if(resampleOccured_){ 
    for(uint i = 0; i < this->nParticles_; i++){
      uint i_prev = this->particleSet_[i]->getParentId();
      if( i_prev != i ){
	landmarkCandidates_[i] = landmarkCandidates_[i_prev];
      }
    }
  }


}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
int FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getGMSize(int i){

  if( i >= 0 && i < this->nParticles_ )
    //return ( maps_[i]->getGaussianCount() );
    return ( this->particleSet_[i]->getData()->getGaussianCount() );
  else
    return -1;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
bool FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
getLandmark(const int i, const int m, 
	    typename TLandmark::Vec &u,
	    typename TLandmark::Mat &S,
	    double &w)
{

  int sz = getGMSize(i);
  if( sz == -1 || (m < 0) || (m >= sz) )
    {
      return false;
    }
  TLandmark *plm;
  //maps_[i]->getGaussian(m, plm, w);
  this->particleSet_[i]->getData()->getGaussian(m, plm, w);
  plm->get(u, S);
  return true;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
setParticlePose(int i, TPose &p){
  
  *(this->particleSet_[i]) = p;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
KalmanFilter* FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getKalmanFilter(){
  return &(kfs_[0]);
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
typename FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::TimingInfo* FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
getTimingInfo(){

  timer_predict_.elapsed(timingInfo_.predict_wall, timingInfo_.predict_cpu);
  timer_mapUpdate_.elapsed(timingInfo_.mapUpdate_wall, timingInfo_.mapUpdate_cpu);
  timer_dataAssoc_.elapsed(timingInfo_.dataAssoc_wall, timingInfo_.dataAssoc_cpu);
  timer_mapUpdate_KF_.elapsed(timingInfo_.mapUpdate_KF_wall, timingInfo_.mapUpdate_KF_cpu);
  timer_weighting_.elapsed(timingInfo_.weighting_wall, timingInfo_.weighting_cpu);
  timer_mapManage_.elapsed(timingInfo_.mapManage_wall, timingInfo_.mapManage_cpu);
  timer_particleResample_.elapsed(timingInfo_.particleResample_wall, timingInfo_.particleResample_cpu);
  
  return &timingInfo_;
}

}

#endif
