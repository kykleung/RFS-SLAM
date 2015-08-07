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

#ifndef RBPHDFILTER_HPP
#define RBPHDFILTER_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/multi_array.hpp>
#include <Eigen/Core>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <list>

#include "GaussianMixture.hpp"
#include "KalmanFilter.hpp"
#include "MurtyAlgorithm.hpp"
#include "ParticleFilter.hpp"
#include "PermutationLexicographic.hpp"
#include "Timer.hpp"
#include "misc/MemProfile.hpp"

namespace rfs{

/**
 *  \class RBPHDFilter
 *  \brief Rao-Blackwellized Probability Hypothesis Density Filter class
 *  
 *  This class implements the Rao-Bloackwellized Probability Hypothesis Density
 *  filter. The constructor of this class will internally instantiate the 
 *  process model for both the robot and landmarks, the measurement model, 
 *  and the Kalman filter. Users have access to these through pointers that 
 *  can be obtained by calling the appropraite get function.
 *
 *  \tparam RobotProcessModel A robot process model derived from ProcessModel
 *  \tparam LmkProcessModel A landmark process model derived from ProcessModel
 *  \tparam MeasurementModel A sensor model derived from MeasurementModel
 *  \tparam KalmanFilter A Kalman filter that uses LmkProcessModel and MeasurementModel
 *  \author Keith Leung, Felipe Inostroza
 */

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
class RBPHDFilter : public ParticleFilter<RobotProcessModel, MeasurementModel,
					  GaussianMixture< typename MeasurementModel::TLandmark > >
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef typename RobotProcessModel::TState TPose; /**< \brief Robot pose */
  typedef typename RobotProcessModel::TInput TInput; /**< \brief Process model input */
  typedef typename MeasurementModel::TLandmark TLandmark; /**< \brief Landmark position */
  typedef typename MeasurementModel::TMeasurement TMeasurement; /**< \brief Sensor measurement */
  typedef GaussianMixture<TLandmark> TGM; /**< \brief Gaussian mixture for landmarks */
  typedef typename TGM::Gaussian TGaussian; /**< \brief Landmark Gaussian  */

  /** 
   * \brief Configurations for this RBPHDFilter 
   */
  struct Config{
    
    /**  New birth Gaussians are set with this weight */
    double birthGaussianWeight_;   

    /**  Number of supporting measurements required to generate a birth Gaussian */
    uint birthGaussianMeasurementCountThreshold_;

    /**  Number of checks before candidate birth Gaussian is removed if there are not enough supporting measurements */
    uint birthGaussianMeasurementCheckThreshold_;

    /**  Mahalanobis distance threshold less than which a measurement counts towards supporting a birth Guassian */
    double birthGaussianMeasurementSupportDist_;

    /**  Number of latest observations below or equal to which birth Gaussians are added regardless or measurement evidence.
     *   This is to accomodate for the case where the robot enters an area sparesly populated with landmarks.
     */
    uint birthGaussianCurrentMeasurementCountThreshold_;

    /**  New Gaussians are only created during map update if the innovation mahalanobis distance 
	 is less than this threshold */
    double newGaussianCreateInnovMDThreshold_;

    /**  number of map states to use for evaluating particle weight
	 0 => empty-set strategy,
	 1 => single-feature strategy,
	 >1 => multi-feature strategy
    */
    int importanceWeightingEvalPointCount_;

    /** The weight above which a Gaussian's mean is considered as a evaluation point for particle importance factor */
    double importanceWeightingEvalPointGuassianWeight_;

    /** The mahalanobis distance threshold used to determine if a possible meaurement-landmark
     *  pairing is significant to worth considering 
     */
    double importanceWeightingMeasurementLikelihoodMDThreshold_;

    /** Gaussian merging Mahalanobis distance threshold */
    double gaussianMergingThreshold_;

    /** Gaussian merging covariance inflation factor */
    double gaussianMergingCovarianceInflationFactor_;

    /** Gaussian pruning weight threshold, below which Gaussians are eliminated from a Gaussian mixture */
    double gaussianPruningThreshold_;

    /** Minimum number of updates betwen resampling of particles*/
    int minUpdatesBeforeResample_;

    /** Minimum numnber of measurements before resampling of particles*/
    int minMeasurementsBeforeResample_;

    /** Use the particle weighting strategty from Single-cluster PHD Filtering by Lee, et. al. */
    bool useClusterProcess_;

  } config;


  /**
   * \brief Elapsed timing information 
   */
  struct TimingInfo {
    long long predict_wall;
    long long predict_cpu;
    long long mapUpdate_wall;
    long long mapUpdate_cpu;
    long long mapUpdate_kf_wall;
    long long mapUpdate_kf_cpu;
    long long particleWeighting_wall;
    long long particleWeighting_cpu;
    long long mapMerge_wall;
    long long mapMerge_cpu;
    long long mapPrune_wall;
    long long mapPrune_cpu;
    long long particleResample_wall;
    long long particleResample_cpu;
  } timingInfo_;

  /**
   * \class BirthGaussianCandidate
   * \brief Birth Gaussian candidate
   */
  class BirthGaussianCandidate : public TLandmark{
  public:
    uint nSupportingMeasurements;
    uint nChecks;
  };

  /** 
   * Constructor 
   * \param n number of particles
   */
  RBPHDFilter(int n);

  /** Destructor */
  ~RBPHDFilter();

  /** 
   * Get the landmark process model
   * \return pointer to the landmark process model
   */
  LmkProcessModel* getLmkProcessModel();

  /**
   * Predict the robot trajectory using the lastest odometry data
   * \param[in] u input 
   * \param[in] dT Timestep timestep, which the motion model may or may not use;
   * \param[in] useModelNoise use the additive noise for the process model
   * \param[in] useInputNoise use the noise fn the input
   * \param[in] birthGaussianCheck check whether birth Gaussians should be created
   */
  void predict( TInput u, TimeStamp const &dT,
		bool useModelNoise = true,
		bool useInputNoise = false,
		bool birthGaussianCheck = true);

  /**
   * Update the map, calculate importance weighting, sample if necessary, and
   * create new birth Gaussians.
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
   * \param[out] w weight
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

  std::vector< std::list<BirthGaussianCandidate> > birthGaussians_;

  typedef boost::multi_array<double, 2> W_Table; /**< \brief weighting table */
  typedef boost::multi_array<TLandmark*, 2> M_Table; /**< \brief landmark table */
  std::vector<W_Table> wTables_; /**< \brief Weighting table used during map update, one for each thread */
  std::vector<M_Table> mTables_; /**< Table for initiating new landmarks during map update, one for each thread */

  int nThreads_; /**< Number of threads for running map update for all particles */

  std::vector<KalmanFilter> kfs_; /**< Kalman filters (one for each thread)*/
  LmkProcessModel *lmkModelPtr_; /**< pointer to landmark process model */
  double threshold_mahalanobisDistance2_mapUpdate_; /**< Mahalanobis distance squared threshold, above which Kalman Filter update will not occur*/

  /** indices of unused measurement for each particle for creating birth Gaussians */
  std::vector< std::vector<unsigned int> > unused_measurements_; 

  /** count of how many landmarks were in the FOV */
  std::vector<uint> nLandmarksInFOV_;

  unsigned int nUpdatesSinceResample_; /**< Number of updates performed since the last resmaple */
  unsigned int nMeasurementsSinceResample_; /**< Number of measurements processed since the last resample */
  bool resampleOccured_;

  Timer timer_predict_; /**< Timer for prediction step */
  Timer timer_mapUpdate_; /**< Timer for map update */
  Timer timer_mapUpdate_kf_; /**<Timer for the Kalman filter part of map update */
  Timer timer_particleWeighting_; /**< Timer for particle weighting */
  Timer timer_mapMerge_; /**< Timer for map merging */ 
  Timer timer_mapPrune_; /**< Timer for map pruning */
  Timer timer_particleResample_; /**<Timer for particle resampling */
  

  /** 
   * Add birth Gaussians for each particle's map using unused_measurements_
   */ 
  void addBirthGaussians();

  /**
   * Update the map with the measurements in measurements_
   * Existing landmarks with probability of detection > 0 will have their Gaussian
   * mixture weight reduced to account for missed detection.
   * For every landmark-measurement pair with probability of detection > 0,
   * a new landmark will be created. 
   * \param[in] particleIdx index of the particle for which the map will be updated
   */
  void updateMap(const uint particleIdx);

  /** 
   * Importance weighting. Overrides the abstract function in ParticleFilter
   * \param[in] idx particle index
   */
  void importanceWeighting(const uint idx);

  /**
   * Random Finite Set measurement likelihood evaluation
   * \brief The current measurements in measurements_ are used to determine the
   * RFS measurement likelihood given a set of landmarks 
   * \param[in] particleIdx particle for which the likelihood is calcuated
   * \param[in] evalPtIdx indices of evaluation points in maps_[particleIdx]
   * \param[in] evalPtPd probability of detection of evaluation point 
   * \return measurement likelihood
   */
  double rfsMeasurementLikelihood( const int particleIdx, 
				   std::vector<unsigned int> &evalPtIdx,
				   std::vector<double> &evalPtPd );

  /**
   * Calculate the sum of all permutations of measurement likelihood from a likelihood
   * table generated from within rfsMeasurementLikelihood
   * \param[in] likelihoodTab likelihood table generated within rfsMeasurementLikelihood
   * \param[in] Z_NoClutter A vector of measurement indices (columns) to consider in the likelihoodTab 
   * \return sum of all permutations from the given likelihood table
   */
  double rfsMeasurementLikelihoodPermutations( std::vector< double* > &likelihoodTab, 
					       std::vector< int > &Z_NoClutter);

  /** Checks the Gaussian mixture maps for all particles for errors
   *  \return true if there are no errors 
   */
  bool checkMapIntegrity();

  /** \brief Ensure that Weighting tables have sufficient capacity
   *  \param[in] threadNum thread number of table
   *  \param[in] nRows desired number of rows
   *  \param[in] nCols desired number of columns
   */
  void checkWeightingTableSize(unsigned int const threadNum, unsigned int nRows, unsigned int nCols);

};

////////// Implementation //////////

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::RBPHDFilter(int n)
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

  for(int i = 0; i < n; i++){
    this->particleSet_[i]->setData( boost::shared_ptr<TGM>( new TGM() ) );
  }

  for(int i = 0; i < n; i++){
    unused_measurements_.push_back( std::vector<unsigned int>() );
  }
  
  config.birthGaussianWeight_ = 0.25; 
  config.birthGaussianMeasurementCountThreshold_ = 1;
  config.birthGaussianMeasurementCheckThreshold_ = 1;
  config.birthGaussianMeasurementSupportDist_ = 1;
  config.birthGaussianCurrentMeasurementCountThreshold_ = 1;
  config.gaussianMergingThreshold_ = 0.5;
  config.gaussianMergingCovarianceInflationFactor_ = 1.5;
  config.gaussianPruningThreshold_ = 0.2;
  config.importanceWeightingEvalPointCount_ = 8;
  config.importanceWeightingMeasurementLikelihoodMDThreshold_ = 3.0;
  config.newGaussianCreateInnovMDThreshold_ = 0.2;
  config.minUpdatesBeforeResample_ = 1;
  config.minMeasurementsBeforeResample_ = 1;
  
  nUpdatesSinceResample_ = 0;
  nMeasurementsSinceResample_ = 0;

  wTables_.reserve(nThreads_);
  mTables_.reserve(nThreads_);
  double const TABLE_INIT_SIZE = 100;
  for(int i = 0 ; i < nThreads_ ; i++){
    wTables_.push_back( W_Table(boost::extents[TABLE_INIT_SIZE][TABLE_INIT_SIZE]) );
    mTables_.push_back( M_Table(boost::extents[TABLE_INIT_SIZE][TABLE_INIT_SIZE]) ); 
  }

  birthGaussians_.resize(n);
  nLandmarksInFOV_.resize(n);
 
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::~RBPHDFilter(){

  for(int i = 0; i < this->nParticles_; i++){
    this->particleSet_[i]->deleteData();
  }
  delete lmkModelPtr_;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
LmkProcessModel* RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getLmkProcessModel(){
  return lmkModelPtr_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::predict( TInput u,
												 TimeStamp const &dT,
												 bool useModelNoise,
												 bool useInputNoise,
												 bool birthGaussianCheck){
 
  timer_predict_.resume();

  // Add birth Gaussians using pose before prediction
  if( birthGaussianCheck ){
    addBirthGaussians();
  }

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
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::update( std::vector<TMeasurement> &Z){

  
  nUpdatesSinceResample_++;

  this->setMeasurements( Z ); // Z gets cleared after this call, measurements now stored in this->measurements_
  if(this->measurements_.size() == 0)
    return;
  nMeasurementsSinceResample_ += this->measurements_.size();

  const unsigned int startIdx = 0;
  const unsigned int stopIdx = this->nParticles_;
  const unsigned int nZ = this->measurements_.size();

  // make sure any setting changes to the Kalman Filter are set for all threads
  if(nThreads_>1){
    for(int j = 1; j < nThreads_; j++){
      kfs_[j]=kfs_[0];
    }
  }

  if(nThreads_ > 1){
    timer_mapUpdate_.resume();
  }
  #pragma omp parallel 
  {
    ////////// Map Update //////////
    if(nThreads_ == 1){
      timer_mapUpdate_.resume();
    }
    #pragma omp for
    for(unsigned int i = startIdx; i < stopIdx; i++){
      updateMap(i);
    }
    if(nThreads_ == 1){
      timer_mapUpdate_.stop();
    }

    ////////// Particle Weighting ////////// 
    if(!config.useClusterProcess_){
      if(nThreads_ == 1){
	timer_particleWeighting_.resume();
      }
      #pragma omp for
      for(unsigned int i = startIdx; i < stopIdx; i++){
	importanceWeighting(i);
      }
      if(nThreads_ == 1){
	timer_particleWeighting_.stop();
      }
    }
      
    //////////// Merge and prune //////////
    if(nThreads_ == 1){
      timer_mapMerge_.resume();
    }
    #pragma omp for
    for(unsigned int i = startIdx; i < stopIdx; i++){
      this->particleSet_[i]->getData()->merge( config.gaussianMergingThreshold_, 
					       config.gaussianMergingCovarianceInflationFactor_);
    }
    if(nThreads_ == 1){
      timer_mapMerge_.stop();
    }
    
    if(nThreads_ == 1){
      timer_mapPrune_.resume();
    }
    #pragma omp for
    for(unsigned int i = startIdx; i < stopIdx; i++){
      this->particleSet_[i]->getData()->prune( config.gaussianPruningThreshold_ );     
    }
    if(nThreads_ == 1){
      timer_mapPrune_.stop();
    }
  }
  if(nThreads_ > 1){
    timer_mapUpdate_.stop();
  }

  //////////// Particle resampling //////////
  timer_particleResample_.resume();
  resampleOccured_ = false;
  if( nUpdatesSinceResample_ >= config.minUpdatesBeforeResample_ && 
      nMeasurementsSinceResample_ >= config.minMeasurementsBeforeResample_){
    resampleOccured_ = this->resample();
  }

  if( resampleOccured_ ){
    nUpdatesSinceResample_ = 0;
    nMeasurementsSinceResample_ = 0;
  }else{
    this->normalizeWeights();
  }
  timer_particleResample_.stop();

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::updateMap(const uint particleIdx){

    const unsigned int i = particleIdx;
    const unsigned int nZ = this->measurements_.size();
    
    int threadnum = 0;
    #ifdef _OPENMP
    threadnum = omp_get_thread_num();
    #endif

    //---------- 1. setup / book-keeping ----------
   
    const unsigned int nM = this->particleSet_[i]->getData()->getGaussianCount();
    unused_measurements_[i].clear();  
    nLandmarksInFOV_[i] = 0;  
    if(nM == 0){ // No existing landmark case -> flag all measurements as unused and go to next particles
      for(int z = 0; z < nZ; z++){
	unused_measurements_[i].push_back( z );
      }
      return; // goto next particle
    }
    double Pd[nM];
    int landmarkCloseToSensingLimit[nM];

    // For cluster process particle weighting

    double w_km_sum = std::numeric_limits<double>::denorm_min();
    double likelihoodProd = 1;
    if(config.useClusterProcess_){
      if(nThreads_ == 1){
	timer_particleWeighting_.resume();
      }
      for(int m = 0; m < nM; m++){
	w_km_sum += this->particleSet_[i]->getData()->getWeight(m);
      }
      timer_particleWeighting_.stop();
    }

    // nM x nZ table for Gaussian weighting
    checkWeightingTableSize(threadnum, nM, nZ);


    //----------  2. Kalman Filter map update ----------

    if(nThreads_ == 1)
      timer_mapUpdate_kf_.resume();
    const TPose *pose = this->particleSet_[i]->getPose();
    threshold_mahalanobisDistance2_mapUpdate_ = config.newGaussianCreateInnovMDThreshold_ * config.newGaussianCreateInnovMDThreshold_;

    std::vector<double> innovationLikelihood(nZ);
    std::vector<double> innovationMahalanobisDist2(nZ);
    std::vector<TLandmark> lmNew(nZ);

    for(unsigned int m = 0; m < nM; m++){

      TLandmark* lm = this->particleSet_[i]->getData()->getGaussian(m);
      bool isCloseToSensingLimit;
      Pd[m] = this->pMeasurementModel_->probabilityOfDetection( *pose, *lm, 
								isCloseToSensingLimit); 

      if(isCloseToSensingLimit){
	landmarkCloseToSensingLimit[m] = 1;
	Pd[m] = 1;
      }else{
	landmarkCloseToSensingLimit[m] = 0;
      }
      double w_km = this->particleSet_[i]->getData()->getWeight(m);
      double Pd_times_w_km = Pd[m] * w_km;

      if(Pd[m] != 0){

	nLandmarksInFOV_[i]++;

	// RUN KF, create new landmark for likely updates but do not add to map_[i] yet
	// because we cannot determine actual weight until the entire weighting table is
	// filled in
	kfs_[threadnum].correct(*pose, this->measurements_, *lm, lmNew, &innovationLikelihood, &innovationMahalanobisDist2);
	
	for(int z = 0; z < nZ; z++){

	  if ( innovationLikelihood[z] == 0 || innovationMahalanobisDist2[z] > threshold_mahalanobisDistance2_mapUpdate_ ){
	    mTables_[threadnum][m][z] = NULL;
	    wTables_[threadnum][m][z] = 0;
	  }else{
	    mTables_[threadnum][m][z] = new TLandmark( lmNew[z] );
	    wTables_[threadnum][m][z] = Pd_times_w_km * innovationLikelihood[z];
	  }	

	} // z forloop end

      }else{ // Pd = 0
	for(int z = 0; z < nZ; z++){
	  mTables_[threadnum][m][z] = NULL;
	  wTables_[threadnum][m][z] = 0;
	}
      }

    } // m forloop end
    
    // Now calculate the weight of each new Gaussian
    for(int z = 0; z < nZ; z++){

      double clutter = this->pMeasurementModel_->clutterIntensity( this->measurements_[z], nZ );
      double sum = clutter;
      for(unsigned int m = 0; m < nM; m++){
	sum += wTables_[threadnum][m][z];
      }

      if(config.useClusterProcess_){
	likelihoodProd *= sum;
      }

      for(unsigned int m = 0; m < nM; m++){
	wTables_[threadnum][m][z] = wTables_[threadnum][m][z] / sum;
      }
    }

    if(config.useClusterProcess_){
      if(nThreads_ == 1){
	timer_particleWeighting_.resume();
      }
	double prev_particle_i_weight = this->particleSet_[i]->getWeight();
	this->particleSet_[i]->setWeight( exp(w_km_sum) * likelihoodProd * prev_particle_i_weight);
      timer_particleWeighting_.stop();
    }

    if(nThreads_ == 1)
      timer_mapUpdate_kf_.stop();

    // ---------- 3. Add new Gaussians to map  ----------
    // New Gaussians will have indices >= nM 
    for(int m = 0; m < nM; m++){
      for(int z = 0; z < nZ; z++){
	if(mTables_[threadnum][m][z] != NULL && wTables_[threadnum][m][z] > 0){
	  this->particleSet_[i]->getData()->addGaussian( mTables_[threadnum][m][z],  
							 wTables_[threadnum][m][z]);  
	  mTables_[threadnum][m][z] = NULL;
	}
      }
    }

    //----------  4. Determine weights for existing Gaussians (missed detection) ----------
    for(int m = 0; m < nM; m++){
      
      double w_km = this->particleSet_[i]->getData()->getWeight(m);
      double w_k = (1 - Pd[m]) * w_km;

      // For landmarks close to sensing limit
      if (landmarkCloseToSensingLimit[m] == 1 && w_km > config.birthGaussianWeight_){
	double weight_sum_m = 0;
	for(int z = 0; z < nZ; z++){
	  weight_sum_m += wTables_[threadnum][m][z];
	}
	double delta_w = Pd[m] * w_km - weight_sum_m;
	if( delta_w > 0 ){
	  w_k += delta_w; // This is just a heuristic that works well
	  if(w_k > 1)
	    w_k = 1;
	}
      }

      this->particleSet_[i]->getData()->setWeight(m, w_k);
    }

    //----------  5. Identify unused measurements for adding birth Gaussians later ----------
    unused_measurements_[i].clear();
    for(int z = 0; z < nZ; z++){
      bool is_measurement_used = false;
      for(int m = 0; m < nM; m++){
	if (wTables_[threadnum][m][z] != 0){
	  is_measurement_used = true;
	  break;
	}
      }
      if (!is_measurement_used)
	unused_measurements_[i].push_back( z );
    }


  

}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::importanceWeighting(const uint idx){

  
    const uint particleIdx = idx;
    TPose x = *(this->particleSet_[particleIdx]);

    // 1. select evaluation points from highest-weighted Gaussians after update, that are within sensor FOV
    const unsigned int nM = this->particleSet_[particleIdx]->getData()->getGaussianCount();
    int nEvalPoints = config.importanceWeightingEvalPointCount_ > nM ? nM : config.importanceWeightingEvalPointCount_ ;
    std::vector<unsigned int> evalPointIdx;
    std::vector<double> evalPointPd;
    evalPointIdx.reserve(nEvalPoints);
    evalPointPd.reserve(nEvalPoints);
    if( nEvalPoints == 0 ){
      this->particleSet_[particleIdx]->setWeight( std::numeric_limits<double>::denorm_min() );
      return;
    }
    this->particleSet_[particleIdx]->getData()->sortByWeight(); // sort by weight so that we can pick off the top nEvalPoints Gaussians
    for(int m = 0; m < nM; m++){
      TLandmark* plm_temp;
      double w, w_prev;
      bool closeToSensingLim;
      this->particleSet_[particleIdx]->getData()->getGaussian(m, plm_temp, w, w_prev);
      if(w < config.importanceWeightingEvalPointGuassianWeight_)
	break;
      double Pd = this->pMeasurementModel_->probabilityOfDetection(x, *plm_temp, closeToSensingLim);
      if( Pd > 0 ){
	evalPointIdx.push_back(m);
	evalPointPd.push_back(Pd);
      }
      if(nEvalPoints != -1 && evalPointIdx.size() >= nEvalPoints)
	break;
    }
    nEvalPoints = evalPointIdx.size();

    // 2. evaluate sum of Gaussian weights
    double gaussianWeightSumBeforeUpdate = 0;
    double gaussianWeightSumAfterUpdate = 0;
    for(int m = 0; m < nM; m++){
      TLandmark* plm_temp;
      double w, w_prev;
      this->particleSet_[particleIdx]->getData()->getGaussian(m, plm_temp, w, w_prev); // for newly created Gaussians, w_prev = 0
      gaussianWeightSumBeforeUpdate += w_prev;
      gaussianWeightSumAfterUpdate += w;
    }
    
    // 3. evaluate intensity function at eval points and take their product
    double intensityProd_beforeUpdate = 1;
    double intensityProd_afterUpdate = 1;
    for(int e = 0; e < nEvalPoints; e++){

      int p = evalPointIdx[e];
      TLandmark* lm_evalPt;
      double w_temp;
      this->particleSet_[particleIdx]->getData()->getGaussian(p, lm_evalPt, w_temp);

      double intensity_at_evalPt_beforeUpdate = std::numeric_limits<double>::denorm_min();
      double intensity_at_evalPt_afterUpdate = std::numeric_limits<double>::denorm_min();

      for(int m = 0; m < nM; m++){
	TLandmark* plm;
	double w, w_prev;
	this->particleSet_[particleIdx]->getData()->getGaussian(m, plm, w, w_prev);
	// New Gaussians from update will have w_prev = 0
	// Old Gaussians (missed-detection) will not have been updated, but weights will have changed
	double likelihood = plm->evalGaussianLikelihood( *lm_evalPt );
	intensity_at_evalPt_beforeUpdate += w_prev * likelihood; // w_prev for newly created Gaussians are 0
	intensity_at_evalPt_afterUpdate += w * likelihood;
      }
      intensityProd_beforeUpdate *= intensity_at_evalPt_beforeUpdate;
      intensityProd_afterUpdate *= intensity_at_evalPt_afterUpdate;
    }

    // 4. calculate measurement likelihood at eval points
    // note that rfsMeasurementLikelihood uses maps_[particleIdx] which is already sorted by weight
    double measurementLikelihood = rfsMeasurementLikelihood( particleIdx, evalPointIdx, evalPointPd );
    //printf("Particle %d measurement likelihood = %f\n", i, measurementLikelihood);

    // 5. calculate overall weight
    double overall_weight = measurementLikelihood * intensityProd_beforeUpdate / intensityProd_afterUpdate *
      exp( gaussianWeightSumAfterUpdate - gaussianWeightSumBeforeUpdate); 
    
    double prev_weight = this->particleSet_[particleIdx]->getWeight();
    this->particleSet_[particleIdx]->setWeight( overall_weight * prev_weight );
    /*std::cout << "P" << idx << ".w " 
	      << std::setw(15) << intensityProd_beforeUpdate
	      << std::setw(15) << intensityProd_afterUpdate
	      << std::setw(15) << measurementLikelihood
	      << std::setw(15) << overall_weight << std::endl;*/

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
double RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
rfsMeasurementLikelihood( const int particleIdx, 
			  std::vector<unsigned int> &evalPtIdx,
			  std::vector<double> &evalPtPd ){

  // eval points are first nEvalPoints elements of maps_[i], which are already ordered by weight; 

  const int i = particleIdx;
  TPose x = *(this->particleSet_[i]);
  const int nM = evalPtIdx.size();
  const int nZ = this->measurements_.size();
  int nL = nM;
  if(nZ > nL)
    nL = nZ;
  TLandmark* evalPt;
  TLandmark evalPt_copy;
  TMeasurement expected_z;
  double const threshold = config.importanceWeightingMeasurementLikelihoodMDThreshold_ *
    config.importanceWeightingMeasurementLikelihoodMDThreshold_;
  double md2; // Mahalanobis distance squared

  // Create and fill in likelihood table (nM x nZ)
  double** L;
  CostMatrixGeneral likelihoodMatrix(L, nM, nZ);

  for(int m = 0; m < nM; m++){
    
    this->particleSet_[i]->getData()->getGaussian( evalPtIdx[m], evalPt ); // get location of m
    evalPt_copy = *evalPt; // so that we don't change the actual data //
    evalPt_copy.setCov(MeasurementModel::TLandmark::Mat::Zero()); //
    this->pMeasurementModel_->measure( x, evalPt_copy, expected_z); // get expected measurement for m
    double Pd = evalPtPd[m]; // get the prob of detection of m

    for(int n = 0; n < nZ; n++){

      // calculate measurement likelihood with detection statistics
      L[m][n] = expected_z.evalGaussianLikelihood( this->measurements_[n], &md2) * Pd; // new line 
      if( md2 > threshold ){
	L[m][n] = 0;
      }
    }
  }

  // Partition the Likelihood Table and turn into a log-likelihood table
  int nP = likelihoodMatrix.partition();
  double l = 1;
  double const BIG_NEG_NUM = -1000; // to represent log(0)
  double clutter[nZ];
  for(int n = 0; n < nZ; n++ ){
    clutter[n] = this->pMeasurementModel_->clutterIntensity( this->measurements_[n], nZ );
  }

  // Go through each partition and determine the likelihood
  for(int p = 0; p < nP; p++){

    double partition_likelihood = 0;

    unsigned int nCols, nRows;
    double** Cp; 
    unsigned int* rowIdx;
    unsigned int* colIdx;    

    bool isZeroPartition = !likelihoodMatrix.getPartitionSize(p, nRows, nCols);
    bool useMurtyAlgorithm = true;
    if(nRows + nCols <= 8 || isZeroPartition)
      useMurtyAlgorithm = false;
  
    isZeroPartition = !likelihoodMatrix.getPartition(p, Cp, nRows, nCols, rowIdx, colIdx, useMurtyAlgorithm); 
   
    if(isZeroPartition){ // all landmarks in this partition are mis-detected. All measurements are outliers
      
      partition_likelihood = 1;
      for(int r = 0; r < nRows; r++){
	partition_likelihood *= evalPtPd[rowIdx[r]];
      }

      for(int c = 0; c < nCols; c++){
	partition_likelihood *= clutter[colIdx[c]];
      }


    }else{
      // turn the matrix into a log likelihood matrix with detection statistics,
      // and fill in the extended part of the partition

      for(int r = 0; r < nRows; r++){
	for(int c = 0; c < nCols; c++){
	  if(Cp[r][c] == 0)
	    Cp[r][c] = BIG_NEG_NUM;
	  else{
	    Cp[r][c] = log(Cp[r][c]);
	    if(Cp[r][c] < BIG_NEG_NUM)
	      Cp[r][c] = BIG_NEG_NUM;
	  }
	}
      }


      if(useMurtyAlgorithm){ // use Murty's algorithm

	// mis-detections
	for(int r = 0; r < nRows; r++){
	  for(int c = nCols; c < nRows + nCols; c++){
	    if(r == c - nCols)
	      Cp[r][c] = log(1 - evalPtPd[rowIdx[r]]); 
	    else
	      Cp[r][c] = BIG_NEG_NUM;
	  }
	}

	// clutter
	for(int r = nRows; r < nRows + nCols; r++){
	  for(int c = 0; c < nCols; c++){
	    if(r - nRows == c)
	      Cp[r][c] = log(clutter[colIdx[c]]); 
	    else
	      Cp[r][c] = BIG_NEG_NUM;
	  }
	}

	// the lower right corner
	for(int r = nRows; r < nRows + nCols; r++){
	  for(int c = nCols; c < nRows + nCols; c++){
	    Cp[r][c] = 0;
	  }
	}

	Murty murtyAlgo(Cp, nRows + nCols);
	Murty::Assignment a;
	partition_likelihood = 0;
	double permutation_log_likelihood = 0;
	murtyAlgo.setRealAssignmentBlock(nRows, nCols);
	for(int k = 0; k < 200; k++){ 
	  int rank = murtyAlgo.findNextBest(a, permutation_log_likelihood);
	  if(rank == -1 || permutation_log_likelihood < BIG_NEG_NUM)
	    break;
	  partition_likelihood += exp(permutation_log_likelihood);
	}  

      }else{ // use lexicographic ordering

	partition_likelihood = 0;
	double permutation_log_likelihood = 0; 

	uint o[nRows + nCols];
 
	PermutationLexicographic pl(nRows, nCols, true);
	unsigned int nPerm = pl.next(o);
	while( nPerm != 0){
	  permutation_log_likelihood = 0; 
	  for(int a = 0; a < nRows; a++){
	    if(o[a] < nCols){ // detection
	      permutation_log_likelihood += Cp[a][o[a]];
	    }else{ // mis-detection
	      permutation_log_likelihood += log(1 - evalPtPd[rowIdx[a]]); 
	    }
	  }
	  for(int a = nRows; a < nRows + nCols; a++){ // outliers
	    if(o[a] < nCols){
	      permutation_log_likelihood += log(clutter[colIdx[o[a]]]);
	    }
	  }
	  partition_likelihood += exp(permutation_log_likelihood);
	  nPerm = pl.next(o);
	}

      } // End lexicographic ordering
    
    } // End non zero partition

    l *= partition_likelihood;

  } // End partitions
 
  return (l / this->pMeasurementModel_->clutterIntensityIntegral( nZ ) );
}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::addBirthGaussians(){

  for(int i = 0; i < this->nParticles_; i++){

    if(resampleOccured_){ 
      uint i_prev = this->particleSet_[i]->getParentId();
      if( i_prev != i ){
	unused_measurements_[i] = unused_measurements_[i_prev];
	birthGaussians_[i] = birthGaussians_[i_prev];
      }
    }
    
    while( unused_measurements_[i].size() > 0){
     
      // get measurement
      int unused_idx = unused_measurements_[i].back();
      TMeasurement unused_z = this->measurements_[unused_idx];
      unused_measurements_[i].pop_back();

      // check to see if measurement corresponds closely with a Landmark in the birth Gaussian list
      bool isNewBirthGaussianCandidate = true;
      
      for( typename std::list<BirthGaussianCandidate>::iterator it = birthGaussians_[i].begin();
	   it != birthGaussians_[i].end(); it++ ){

	TMeasurement z_exp;
	TPose x = *(this->particleSet_[i]); 
	this->pMeasurementModel_->measure(x, *it, z_exp);
	double d2 = z_exp.mahalanobisDist2( unused_z );
	if(d2 <= config.birthGaussianMeasurementSupportDist_ * config.birthGaussianMeasurementSupportDist_){
	  kfs_[0].correct(x, unused_z, *it, *it);
	  (it->nSupportingMeasurements)++;
	  isNewBirthGaussianCandidate = false;
	  break;
	}
      }
	
      if(isNewBirthGaussianCandidate){

	// use inverse measurement model to get landmark
	TPose robot_pose;
	BirthGaussianCandidate c;
	c.nSupportingMeasurements = 1;
	c.nChecks = 0;
	// this->particleSet_[i]->getPose(robot_pose);
	robot_pose = *(this->particleSet_[i]); 
	this->pMeasurementModel_->inverseMeasure( robot_pose, unused_z, c );

	if(config.birthGaussianMeasurementCountThreshold_ == 1 ||
	   nLandmarksInFOV_[i] <= config.birthGaussianCurrentMeasurementCountThreshold_){
	  // add birth landmark to Gaussian mixture (last param = true to allocate mem)
	  this->particleSet_[i]->getData()->addGaussian( &c, config.birthGaussianWeight_, true);	  
	}else{
	  birthGaussians_[i].push_back(c);
	}
	
      }  

    } // iterate unused measurements end

    // Check through candidate list to see if any candidates should be made into a real birth Gaussian
    for( typename std::list<BirthGaussianCandidate>::iterator it = birthGaussians_[i].begin();
	 it != birthGaussians_[i].end(); it++ ){
      it->nChecks++;
      while( it->nSupportingMeasurements >= config.birthGaussianMeasurementCountThreshold_ || 
	     it->nChecks > config.birthGaussianMeasurementCheckThreshold_ ||
	     nLandmarksInFOV_[i] <= config.birthGaussianCurrentMeasurementCountThreshold_){
	if(it->nSupportingMeasurements >= config.birthGaussianMeasurementCountThreshold_){
	  this->particleSet_[i]->getData()->addGaussian( &(*it), config.birthGaussianWeight_, true);
	}else if( nLandmarksInFOV_[i] <= config.birthGaussianCurrentMeasurementCountThreshold_){
	  this->particleSet_[i]->getData()->addGaussian( &(*it), config.birthGaussianWeight_, true);
	}
	it = birthGaussians_[i].erase( it );
	if( it != birthGaussians_[i].end() ){
	  it->nChecks++;
	}else{
	  break;
	}
      }
    }
    
  } // iterate particles end

}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
bool RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::checkMapIntegrity(){

  for( int i = 0; i < this->nParticles_; i++ ){

    unsigned int nM = this->particleSet_[i]->getData()->getGaussianCount();
    for(int m = 0; m < nM; m++){

      TLandmark *lm;
      typename TLandmark::Vec lm_x;
      typename TLandmark::Mat lm_S;
      double w;

      this->particleSet_[i]->getData()->getGaussian(m, lm, w);
      if( lm != NULL ){
	lm->get(lm_x, lm_S);

	bool vecError = false;
	for( int c = 0; c < lm_x.rows(); c++){
	  if(lm_x(c) != lm_x(c)){
	    vecError = true;
	    break;
	  }
	}
	if(vecError){
	  printf("particle %d, landmark index %d, vector error\n", i, m);
	  std::cout << lm_x << std::endl;
	  return false;
	}

	bool matError = false;
	for( int r = 0; r < lm_S.rows(); r++){
	  for( int c = 0; c < lm_S.cols(); c++){
	    if(lm_S(r,c) != lm_S(r,c)){
	      matError = true;
	      break;
	    }
	  }
	}
	lm_x.setOnes();
	double posDefCheck = lm_x.transpose() * lm_S * lm_x;
	if( posDefCheck != posDefCheck || posDefCheck <= 0){
	  matError = true;
	}
	if(matError){
	  printf("particle %d, landmark index %d, covariance error\n", i, m);
	  std::cout << lm_S << std::endl;
	  return false;
	}

	if( w != w ){
	  printf("particle %d, landmark index %d, w = %f\n", i, m, w);
	  return false;
	}

      }

    }

  }

  return true;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
int RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getGMSize(int i){

  if( i >= 0 && i < this->nParticles_ )
    return ( this->particleSet_[i]->getData()->getGaussianCount() );
  else
    return -1;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
bool RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
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
    this->particleSet_[i]->getData()->getGaussian(m, plm, w);
    plm->get(u, S);
    return true;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
setParticlePose(int i, TPose &p){
  
  *(this->particleSet_[i]) = p;

}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
KalmanFilter* RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getKalmanFilter(){
  return &(kfs_[0]);
}



template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::checkWeightingTableSize( unsigned int const threadNum, unsigned int nRows, unsigned int nCols){

  bool resize = false;
  if( wTables_[threadNum].shape()[0] <= nRows ){
    nRows *= 1.2; // expand number of rows by a factor of 1.2
    resize = true;
  }else{
    nRows = wTables_[threadNum].shape()[0];
  }
  if( wTables_[threadNum].shape()[1] <= nCols ){
    nCols *= 1.2; // expand number of cols by a factor of 1.2
    resize = true;
  }else{
    nCols = wTables_[threadNum].shape()[1];
  }
  if(resize){
    wTables_[threadNum].resize( boost::extents[nRows][nCols] );
    mTables_[threadNum].resize( boost::extents[nRows][nCols] );
  }
  
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
typename RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::TimingInfo* RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
  getTimingInfo(){

  timer_predict_.elapsed(timingInfo_.predict_wall, timingInfo_.predict_cpu);
  timer_mapUpdate_.elapsed(timingInfo_.mapUpdate_wall, timingInfo_.mapUpdate_cpu);
  timer_mapUpdate_kf_.elapsed(timingInfo_.mapUpdate_kf_wall, timingInfo_.mapUpdate_kf_cpu);
  timer_particleWeighting_.elapsed(timingInfo_.particleWeighting_wall, timingInfo_.particleWeighting_cpu);
  timer_mapMerge_.elapsed(timingInfo_.mapMerge_wall, timingInfo_.mapMerge_cpu);
  timer_mapPrune_.elapsed(timingInfo_.mapPrune_wall, timingInfo_.mapPrune_cpu);
  timer_particleResample_.elapsed(timingInfo_.particleResample_wall, timingInfo_.particleResample_cpu);
  
  return &timingInfo_;
}

}


#endif
