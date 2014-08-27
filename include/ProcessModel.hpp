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

#ifndef PROCESSMODEL_HPP
#define PROCESSMODEL_HPP

#include "Pose.hpp"
#include "Measurement.hpp"

namespace rfs{

/**
 * \class ProcessModel
 * An abstract class for defining process models
 * \f[ \mathbf{x}_k = \mathbf{g}(\mathbf{x}_{k-1}, \mathbf{u}_k) + \boldsymbol{\delta}, \quad \delta \sim (\mathbf{0}, \mathbf{Q}) \f]
 * where, \f$\mathbf{x}_k\f$ is the updated state,
 * \f$ \mathbf{x}_{k-1} \f$ is the previous state,
 * \f$ \mathbf{u}_k \f$ is the input,
 * \f$ \boldsymbol{\delta} \f$ is the additive white process noise, 
 * \f$ \mathbf{Q} \f$ is the process noise covariance
 * \brief An abstract class for defining process models
 * \tparam StateType RandomVector derived type for state \f$\mathbf{x}\f$ 
 * \tparam InputType Measurement derived for process input \f$\mathbf{u}\f$
 * \author Keith Leung
 */
template<class StateType, class InputType>
class ProcessModel
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef StateType TState;
  typedef InputType TInput;

  /** Default constructor */
  ProcessModel():
    inputNoiseDefined_(false)
  {}

  /** 
   * Constructor 
   * \param[in] Q covariance for the additive zero mean white Gaussian noise for this model 
   */
  ProcessModel(typename StateType::Mat &Q){
    setNoise(Q);
  }

  /** Default destructor */
  ~ProcessModel(){};

  /** 
   * Set the additive zero mean Gaussian noise covariance matrix 
   * \param[in] Q noise covariance
   */
  void setNoise(typename StateType::Mat &Q){
    Q_ = Q;
    inputNoiseDefined_ = true;
  }

  /** 
   * Get the additive zero mean Gaussian noise covariance matrix 
   * \param[in] Q noise covariance
   */
  void getNoise(typename StateType::Mat &Q){
    if(inputNoiseDefined_)
      Q = Q_;
  }

  /** 
   * Abstract function for determining pose at time-step k from pose at time-step k-1
   * This must be implemented in a derived class
   * \param[out] s_k \f$\mathbf{x}_k\f$ pose at current time-step k. This can be the same object as s_km (to update in place)
   * \param[in] s_km \f$\mathbf{x}_{k-1}\f$ pose at previous time-step k-1
   * \param[in] input_k \f$\mathbf{u}_k\f$ input to the process model
   * \param[in] dT size of time-step
   */
  virtual void step( StateType &s_k, StateType &s_km, 
		     InputType &input_k , TimeStamp const &dT) = 0;		     

  /**
   * Sample the process noise to predict the pose at k from k-1
   * \note This function can be overwritten in derived classes for implementing
   * other user-defined sampling methods. 
   * \warning This function does not check that the noise covariance matrices
   * are valid (i.e., semi-positive definite)
   * \param[out] s_k \f$\mathbf{x}_k\f$ sampled pose at current time-step k. This can be the same object as s_km (to update in place). 
   * \param[in] s_km \f$\mathbf{x}_{k-1}\f$ pose at previous time-step k-1
   * \param[in] input_k \f$\mathbf{u}_k\f$ input to process model. If using useInputWhiteGaussianNoise, the assoicated noise needs 
   * to be manually set according to dT.
   * \param[in] dT size of time-step
   * \param[in] useAdditiveWhiteGaussianNoise if true, the output includes 
   * the zero-mean additive white Gaussian noise specified for this ProcessModel
   * \param[in] useInputWhiteGaussianNoise if true, the output includes
   * the noise specified in the input vector, and assumes that it is zero-mean white
   * Gaussian noise.
   */
  virtual void sample( StateType &s_k, StateType &s_km, 
		       InputType &input_k, TimeStamp const &dT,
		       bool useAdditiveWhiteGaussianNoise = true,
		       bool useInputWhiteGaussianNoise = false ){
    
    if(useInputWhiteGaussianNoise){

      InputType in;
      input_k.sample(in); // noise of in needs to be defined according to dT outside this function
      step( s_k, s_km, in, dT );

    }else{
    
      step( s_k, s_km, input_k, dT );

    }
    
    if( useAdditiveWhiteGaussianNoise && Q_ != StateType::Mat::Zero() ){

      s_k.setCov(Q_);
      s_k.sample();
    }
  }

protected:
  
  /** Covariance matrix for zero mean white Gaussian noise */
  typename StateType::Mat Q_;

  /** Flag to indicate if Q_ has been assigned a value */
  bool inputNoiseDefined_;

};

/**
 * \class StaticProcessModel
 * A template process model class with not inputs, used for landmarks
 * \brief A template process model class with not inputs, used for landmarks
 * \author Keith Leung
 */
template< class StateType >
class StaticProcessModel : public ProcessModel< StateType, NullInput>
{

public:

  /** Default constructor */
  StaticProcessModel(){}

  /** Constructor
   *  \param Q additive zero-mean Gaussian noise for this model
   */ 
  StaticProcessModel(typename StateType::Mat &Q): 
    ProcessModel< StateType, NullInput>(Q){
  }

  /** Default destructor */
  ~StaticProcessModel(){}

  /** 
   * Define the step function required by a derived process model and
   * determine the pose at time-step k from pose at time-step k-1
   * \param[out] s_k pose at current time-step k
   * \param[in] s_km pose at previous time-step k-1
   * \param[in] input_k input to process model
   * \param[in] dT size of time-step
   */
  void step( StateType &s_k, StateType &s_km, 
	     NullInput &input_k , TimeStamp const &dT){
    
    if( this->inputNoiseDefined_ ){
      typename StateType::Vec x;
      typename StateType::Mat S;
      s_km.get(x, S, t_);
      S += (this->Q_);
      t_ += dT;
      s_k.set(x, S, t_);
    }else{
      s_k = s_km;
    }
  }
  
  /** 
   * Step function to allow for no inputs to the process model
   * \param[out] s_k State at current time-step k.
   * \param[in] s_km State at previous time-step k-1
   * \param[in] dT size of time-step
   */		     
  void staticStep( StateType &s_k, StateType &s_km, TimeStamp const &dT){
    NullInput input;
    step(s_k , s_km , input , dT);
  }

private:

  TimeStamp t_; 

};

}

#endif
