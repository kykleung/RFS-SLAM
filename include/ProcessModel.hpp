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

// Classes related to motion model
// Keith Leung 2013

#ifndef PROCESSMODEL_HPP
#define PROCESSMODEL_HPP

#include "Measurement.hpp"
#include "Pose.hpp"

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
   * Abstract function for determining pose at time-step k from pose at time-step k-1
   * This must be implemented in a derived class
   * \param[out] s_k \f$\mathbf{x}_k\f$ pose at current time-step k. This can be the same object as s_km (to update in place)
   * \param[in] s_km \f$\mathbf{x}_{k-1}\f$ pose at previous time-step k-1
   * \param[in] input_k \f$\mathbf{u}_k\f$ input to the process model
   * \param[in] dT size of time-step
   */
  virtual void step( StateType &s_k, StateType &s_km, 
		     InputType &input_k , double const dT = 0 ) = 0;		     

  /**
   * Sample the process noise to predict the pose at k from k-1
   * \note This function can be overwritten in derived classes for implementing
   * other user-defined sampling methods. 
   * \warning This function does not check that the noise covariance matrices
   * are valid (i.e., semi-positive definite)
   * \param[out] s_k \f$\mathbf{x}_k\f$ sampled pose at current time-step k. This can be the same object as s_km (to update in place). 
   * \param[in] s_km \f$\mathbf{x}_{k-1}\f$ pose at previous time-step k-1
   * \param[in] input_k \f$\mathbf{u}_k\f$ input to process model
   * \param[in] dT size of time-step
   * \param[in] useAdditiveWhiteGaussianNoise if true, the output includes 
   * the zero-mean additive white Gaussian noise specified for this ProcessModel
   * \param[in] useInputWhiteGaussianNoise if true, the output includes
   * the noise specified in the input vector, and assumes that it is zero-mean white
   * Gaussian noise.
   */
  virtual void sample( StateType &s_k, StateType &s_km, 
		       InputType &input_k, double const dT = 0,
		       bool useAdditiveWhiteGaussianNoise = true,
		       bool useInputWhiteGaussianNoise = false ){
    
    if(useInputWhiteGaussianNoise){

      InputType in;
      input_k.sample(in);
      //RandomVecMathTools<InputType>::sample(input_k, in);

      step( s_k, s_km, in, dT );

    }else{
    
      step( s_k, s_km, input_k, dT );

    }
    
    if( useAdditiveWhiteGaussianNoise && Q_ != StateType::Mat::Zero() ){

      //typename StateType::Vec x_k;
      //double t;
      s_k.setCov(Q_);
      //RandomVecMathTools<StateType>::sample(x_k, Q_, L_, t, s_k);
      s_k.sample();
    }
  }

protected:
  
  /** Covariance matrix for zero mean white Gaussian noise */
  typename StateType::Mat Q_;

  /** Lower triangular part of Cholesky decomposition on S_zmgn */
  //typename StateType::Mat L_;

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
	     NullInput &input_k , double const dT = 1 ){
    
    if( this->inputNoiseDefined_ ){
      typename StateType::Vec x;
      typename StateType::Mat S;
      double t;
      s_km.get(x, S, t);
      S += (this->Q_ * dT * dT);
      t += dT;
      s_k.set(x, S, t);
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
  void staticStep( StateType &s_k, StateType &s_km, double const dT = 1){
    NullInput input;
    step(s_k , s_km , input , dT);
  }

};


/************* 2d Odometry Motion Model *************/

/**
 * \class OdometryMotionModel2d
 * A 2d odometry motion model with translational and rotational
 *        displacement
 * The 2d motion model is as follows: 
 * \f[ \mathbf{x}_k = \begin{bmatrix} \mathbf{p} \\ \theta \end{bmatrix}_k =
 * \begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_k =
 * \mathbf{g}(\mathbf{x}_{k-1}, \mathbf{u}_k) + \boldsymbol{\delta} = 
 * \mathbf{g}\left(\begin{bmatrix} \mathbf{p} \\ \theta \end{bmatrix}_{k-1}, 
 * \begin{bmatrix} \delta \mathbf{p} \\ \delta\theta \end{bmatrix}_k \right) + \boldsymbol{\delta} = 
 * \mathbf{g}\left(\begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_{k-1}, 
 * \begin{bmatrix} \delta x \\ \delta y \\ \delta\theta \end{bmatrix}_k \right) + \boldsymbol{\delta},
 * \boldsymbol\delta \sim (\mathbf{0}, \mathbf{Q})
 * \f]
 * Using rotation matrices to represent the rotation, the function \f$\mathbf{g}\f$
 * composes of the following:
 * \f[
 * \mathbf{p}_k = \mathbf{p}_{k-1} + \mathbf{C}_{k-1}^{\mathrm{T}}(\theta_{k-1}) \delta \mathbf{p}_{k-1}
 * \f]
 * \f[
 * \mathbf{C}_k(\mathbf{\theta_k}) = \mathbf{C}_{k,k-1}(\mathbf{\delta\theta_k}) \mathbf{C}_{k-1}(\mathbf{\theta_{k-1}})
 * \f]
 * \brief A 2d odometry motion model with translational and rotational
 *        displacement
 * \note Currently the updated state from step does not contain valid
 * covariance information because it is not needed by the RBPHDFilter
 * \author Keith Leung
 */
class OdometryMotionModel2d : public ProcessModel< Pose2d, Odometry2d >
{
public:

  /** Default constructor */
  OdometryMotionModel2d();

  /** Constructor with process noise input 
   * \param S additive zero-mean white Gaussian noise covariance matrix
   */
  OdometryMotionModel2d( Pose2d::Mat &Q );

  /** Default destructor */
  ~OdometryMotionModel2d();

  /** 
   * This overrides the virtual function in the parent class for
   * determining the pose at time-step k from pose at time-step k-1.
   * The 2d motion model is as follows: 
   * \f[ \mathbf{x}_k = \begin{bmatrix} \mathbf{p} \\ \theta \end{bmatrix}_k =
   * \begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_k =
   * \mathbf{g}(\mathbf{x}_{k-1}, \mathbf{u}_k) + \boldsymbol{\delta} = 
   * \mathbf{g}\left(\begin{bmatrix} \mathbf{p} \\ \theta \end{bmatrix}_{k-1}, 
   * \begin{bmatrix} \delta \mathbf{p} \\ \delta\theta \end{bmatrix}_k \right) + \boldsymbol{\delta} = 
   * \mathbf{g}\left(\begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_{k-1}, 
   * \begin{bmatrix} \delta x \\ \delta y \\ \delta\theta \end{bmatrix}_k \right) + \boldsymbol{\delta},
   * \boldsymbol\delta \sim (\mathbf{0}, \mathbf{Q})
   * \f]
   * Using rotation matrices to represent the rotation, the function \f$\mathbf{g}\f$
   * composes of the following:
   * \f[
   * \mathbf{p}_k = \mathbf{p}_{k-1} + \mathbf{C}_{k-1}^{\mathrm{T}}(\theta_{k-1}) \delta \mathbf{p}_{k-1}
   * \f]
   * \f[
   * \mathbf{C}_k(\mathbf{\theta_k}) = \mathbf{C}_{k,k-1}(\mathbf{\delta\theta_k}) \mathbf{C}_{k-1}(\mathbf{\theta_{k-1}})
   * \f]
   * \note Currently the updated state from step does not contain valid
   * covariance information because it is not needed by the RBPHDFilter
   * \param[out] s_k pose at current time-step k
   * \param[in] s_km pose at previous time-step k-1
   * \param[in] input_k input to process model
   * \param[in] dT size of time-step (not used)
   */
  void step( Pose2d &s_k, Pose2d &s_km, Odometry2d &input_k, 
	     double const dT=0);
  
};



/**
 * \class OdometryMotionModel1d
 * A 1d odometry motion model with translational displacement
 * The 1d model is as follows:
 * \f[ x_{k} = x_{k-1} + \delta x_k\f]
 * \brief A 1d odometry motion model with translational displacement
 * \note Currently the updated state from step does not contain valid
 * covariance information because it is not needed by the RBPHDFilter
 * \author Keith Leung
 */
class OdometryMotionModel1d : public ProcessModel< Pose1d, Odometry1d >
{
public:

  /** Default constructor */
  OdometryMotionModel1d(){}

  /** Constructor with process noise input 
   * \param[in] Q additive zero-mean white Gaussian noise variance
   */
  OdometryMotionModel1d( Pose1d::Mat &Q );

  /** Default destructor */
  ~OdometryMotionModel1d(){}

   /** 
   * This overrides the virtual function in the parent class for
   * determining the position at time-step k from position at time-step k-1
   * The 1d model is as follows:
   * \f[ x_{k} = x_{k-1} + \delta x_k\f]
   * \note Currently the updated state from step does not contain valid
   * covariance information because it is not needed by the RBPHDFilter
   * \param[out] s_k position at current time-step k
   * \param[in] s_km position at previous time-step k-1
   * \param[in] input_k input to process model
   * \param[in] dT size of time-step (not used)
   */
  void step( Pose1d &s_k, Pose1d &s_km, Odometry1d &input_k, 
	     double const dT=0);

};


#endif
