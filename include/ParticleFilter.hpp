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

#ifndef PARTICLE_FILTER_HPP
#define PARTICLE_FILTER_HPP

#include "Particle.hpp"
#include "ProcessModel.hpp"
#include "MeasurementModel.hpp"
#include <vector>

namespace rfs{

/** 
 * \class ParticleFilter
 * \brief A class containing functions for implementing the particle filter
 * \tparam ProcessModel class for the process model for propagating particles
 * \tparam MeasurementModel class for the measurement model for updateing particles
 * \author Keith Leung
 */
template< class ProcessModel, class MeasurementModel, class ParticleExtraData = int>
class ParticleFilter
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef typename ProcessModel::TState TPose;
  typedef typename ProcessModel::TInput TInput;
  typedef typename MeasurementModel::TMeasurement TMeasure;
  typedef Particle<TPose, ParticleExtraData>* pParticle;
  typedef std::vector<pParticle> TParticleSet;


  /** 
   * \brief Configurations for this ParticleFilter
   */
  struct Config{
    
  }PFconfig;


  /** Defailt constructor */
  ParticleFilter();

  /** 
   * Constructor 
   * \param[in] n number of particles
   * \param[in] initState if not NULL, all particles will take this initial state
   */ 
  ParticleFilter(int n, TPose* initState = NULL);

  /** Default destructor */
  ~ParticleFilter();

  /**
   * Add particles to the particle set
   * \param[in] n the number of particles to add
   * \param[in] initState the initial state of the particles to copy to the new particles
   * \param[in] initWeight the initial importance weighting to assign to the new particles
   * \return the number of particles currently used by the filter after adding
   */
  unsigned int addParticles(int n, TPose* initState = NULL, double initWeight = 1);

  /**
   * Make copies of a particle
   * \param[in] idx the index of the particle to copy
   * \param[in] n the number of copies to make
   * \param[in] weight if non-negative, this weight is assigned to the origin and newly copied particles
   * \return the number of particles currently used by the filter after adding
   */
  unsigned int copyParticle(int idx, int n, double weight = -1);

  /** 
   * Get the process model pointer
   * \return pointer to process model
   */
  ProcessModel* getProcessModel();

  /** 
   * Get the measurement model pointer
   * \return pointer to measurement model
   */
  MeasurementModel* getMeasurementModel();

  /** 
   * Set the measurements for use in importance weight calculations
   * \note The input vector gets cleared
   * \param[in] Z vector container of measurements
   */
  void setMeasurements(std::vector<TMeasure> &Z);

  /** 
   * Propagate particles using the process model
   * \param[in] input to the process model
   * \param[in] dT time-step of input 
   * \param[in] useModelNoise use the additive noise for the process model
   * \param[in] useInputNoise use the noise fn the input
   */
  void propagate( TInput &input, TimeStamp const &dT, 
		  bool useModelNoise = true,
		  bool useInputNoise = false);

  /**
   * Calculate and update importance weights for all particles;
   * Derived class needs to implement this method
   */
  virtual void importanceWeighting();

  /** 
   * Get the number of particles
   * \return count
   */
  int getParticleCount();

  /** 
   * Get the pointer to the particle container
   * \return a pointer
   */
  TParticleSet* getParticleSet(){return &particleSet_;}

  /** 
   * Set the effective particle count below which we initiate resampling
   * \param[in] t threshold
   */
  void setEffectiveParticleCountThreshold(double t);

  /**
   * Get the effective particle count threshold
   * \return threshold
   */ 
  double getEffectiveParticleCountThreshold();

  /**
   * Particle resampling using a low variance sampling. 
   * Sampling will not occur if number of effective particles is above effNParticles_t_.
   * \param[in] n number of particles in the resampled set. This must be less than or equal to the nu,ber of particles.
   *        The default value of 0 will keep the same number of particles.
   * \param[in] forceResample Force resampling to occur irrespective of the effective number of particles
   * \return true if resampling occured
   */
  bool resample( unsigned int n = 0, bool forceResample = false );

 

protected:

  int nParticles_; /**< number of particles */
  TParticleSet particleSet_; /**< container for particle pointers */

  ProcessModel* pProcessModel_; /**< Process model pointer */
  MeasurementModel* pMeasurementModel_; /**< Measurement model pointer */
  
  double effNParticles_t_; /**< Effective particle count threshold for resampling */
  double effNParticles_t_percent_; /**< Effective particle count percentage threshold for resampling */

  std::vector<TMeasure> measurements_; /**< Container for measurements to use for update of particle weight */

  /** 
   * Normalize particle weights so that they sum to 1
   */
  void normalizeWeights();
  
};

////////// Implementation //////////

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::
ParticleFilter(){
  nParticles_ = 0;
  pProcessModel_ =  new ProcessModel;
  pMeasurementModel_ =  new MeasurementModel;
}


template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::
ParticleFilter(int n, TPose* initState){
  
  // initiate particles
  nParticles_ = 0;
  addParticles(n, initState, 1);

  pProcessModel_ =  new ProcessModel;
  pMeasurementModel_ =  new MeasurementModel;

  effNParticles_t_ = double(nParticles_)/4.0; // default is 1/4 of n

  // set random seed for particle resampling
  srand48((unsigned int)time(NULL));

}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
unsigned int ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::
addParticles(int n, TPose* initState, double initWeight){
  
  // initiate particles
  nParticles_ += n;
  particleSet_.reserve(nParticles_);
 
  TimeStamp t0;

  bool noInitState = true; 
  if(initState == NULL){
    typename TPose::Vec x0;
    x0.setZero();  
    initState = new TPose(x0, t0);
  }else{
    noInitState = false;
  }

  if(initWeight < 0)
    initWeight = 0;
  for( int i = 0 ; i < n ; i++ ){
    particleSet_.push_back( new Particle<TPose, ParticleExtraData>(i, *initState, initWeight) );
  }

  if( noInitState ){
    delete initState;
  }

  return nParticles_;
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
unsigned int ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::
copyParticle(int idx, int n, double weight){
  
  // initiate particles
  nParticles_ += n;
  particleSet_.reserve(nParticles_);
  if( weight < 0 )
    weight = particleSet_[idx]->getWeight();
  else
    particleSet_[idx]->setWeight(weight);
  for( int i = 0 ; i < n ; i++ ){			
    particleSet_.push_back( new Particle<TPose, ParticleExtraData>() );	       
    int idxNew = particleSet_.size() - 1;
    particleSet_[idxNew]->setWeight( weight );
    particleSet_[idxNew]->setId( idxNew );
    particleSet_[idxNew]->setParentId( particleSet_[idx]->getParentId() );
    particleSet_[idx]->copyStateTo( particleSet_[idxNew] );
    particleSet_[idx]->copyDataTo( particleSet_[idxNew] );
  } 
  return nParticles_;
}


template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::~ParticleFilter(){
  
  for( int i = 0 ; i < nParticles_ ; i++ ){
    delete particleSet_[i];
  }
  delete pProcessModel_;
  delete pMeasurementModel_;
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
ProcessModel* ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::
getProcessModel(){
  return pProcessModel_; 
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
MeasurementModel* ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::
getMeasurementModel(){
  return pMeasurementModel_; 
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
void ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::setMeasurements(std::vector<TMeasure> &Z){
  measurements_.swap(Z);
  Z.clear();
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
void ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::propagate( TInput &input, 
										   TimeStamp const &dT,
										   bool useModelNoise,
										   bool useInputNoise){
  TPose x_km, x_k;
  for( int i = 0 ; i < nParticles_ ; i++ ){
    particleSet_[i]->getPose( x_km );
    pProcessModel_->sample( x_k, x_km, input, dT, useModelNoise, useInputNoise);
    particleSet_[i]->setPose( x_k );
  } 
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
void ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::importanceWeighting(){
  return;
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
void ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::normalizeWeights(){
  
  double sum = 0;
  for( int i = 0; i < nParticles_; i++ ){
    sum += particleSet_[i]->getWeight();
  }
  for( int i = 0; i < nParticles_; i++ ){
    particleSet_[i]->setWeight( particleSet_[i]->getWeight() / sum );
  }

}


template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
int ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::getParticleCount(){
  return nParticles_;
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
void ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::
setEffectiveParticleCountThreshold(double t){
  effNParticles_t_ = t;
  effNParticles_t_percent_ = t / nParticles_;
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
double ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::
getEffectiveParticleCountThreshold(){
  return effNParticles_t_;
}

template< class ProcessModel, class MeasurementModel, class ParticleExtraData>
bool ParticleFilter<ProcessModel, MeasurementModel, ParticleExtraData>::resample( unsigned int n, bool forceResample ){

  normalizeWeights(); // sum of all particle weights is now 1

  // Check effective number of particles
  if(!forceResample){
    double sum_of_weight_squared = 0;
    for( int i = 0; i < nParticles_; i++ ){
      double w_i = particleSet_[i]->getWeight();
      sum_of_weight_squared += (w_i * w_i); // and divide by 1
    }
    double nEffParticles_ = 1.0 / sum_of_weight_squared;
    if( nEffParticles_ > effNParticles_t_ && nEffParticles_ / nParticles_ > effNParticles_t_percent_){
      return false; // no resampling
    }
  }

  if( n == 0 || n > nParticles_ )
    n = nParticles_; // number of particles to sample

  // Sampler settings
  double randomNum_0_to_1 = drand48();
  unsigned int idx = 0;
  const double sample_interval = 1.0 / double(n); 
  const double sampler_offset = sample_interval * randomNum_0_to_1;
  double sample_point = sampler_offset;
  double cumulative_weight = particleSet_[idx]->getWeight();

  // book-keeping
  std::vector<char> flag_particle_sampled (nParticles_, 0); 
  std::vector<unsigned int> sampled_idx (n, 0);
  
  // Sample
  for( int i = 0; i < n; i++ ){

    while( sample_point > cumulative_weight ){
      // particle[idx] not sampled
      idx++;
      cumulative_weight += particleSet_[idx]->getWeight();
    }
    // particle[idx] sampled
    sampled_idx[i] = idx;
    flag_particle_sampled[idx] = 1;
    sample_point += sample_interval;
  }

  // Do the actual data copying
  unsigned int idx_prev = 0;
  unsigned int next_unsampled_idx = 0;
  for( int i = 0; i < n; i++ ){
    
    bool firstTime = true;
    idx = sampled_idx[i]; // particle[idx] was sampled 
    
    if(i > 0 && idx == idx_prev ){
      firstTime = false;
    }
    idx_prev = idx;

    // cases:
    // 1. idx < n AND idx appears for first time -> do nothing
    // 2. idx < n AND it is not the first time that idx appears -> make copy
    // 3. idx > n AND idx appears for first time -> make copy
    // 4. idx > n AND it is not first time that idx appears -> make copy

    if( idx < n && firstTime){ // case 1
      particleSet_[idx]->setParentId( particleSet_[idx]->getId() );
    }else{ // case 2, 3, 4

      if( next_unsampled_idx < nParticles_ ){
	while( flag_particle_sampled[next_unsampled_idx] == 1)
	  next_unsampled_idx++;
      }

      particleSet_[next_unsampled_idx]->setParentId( particleSet_[idx]->getId() );
      particleSet_[idx]->copyStateTo( particleSet_[next_unsampled_idx] );
      particleSet_[next_unsampled_idx]->deleteData();
      particleSet_[idx]->copyDataTo( particleSet_[next_unsampled_idx] );

      next_unsampled_idx++;
    }      
  }

  // Delete all particles with idx >= n
  for( int i = n; i < nParticles_; i++ ){
    delete particleSet_[i];
  }
  nParticles_ = n;
  particleSet_.resize(nParticles_);

  // Reset weight of all particles
  for( int i = 0; i < n; i++ ){
    particleSet_[i]->setWeight(1);
  }
  
  return true;

}

}

#endif
