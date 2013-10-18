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

// Particle Class used in the RB-PHD-Filter
// Keith Leung 2013

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

/** 
 *  \class Particle
 *  \brief A class for a particle for the particle filter
 *  \tparam PoseType RandomVec derived class to represent the state of a particle
 *  \author Keith Leung
 */
template< class PoseType >
class Particle
{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef PoseType tPose;

  /** Default constructor */
  Particle();

  /** Constructor 
   * \param[in] id particle id
   * \param[in] x_k_i particle pose
   * \param[in] w particle weight
   */
  Particle( unsigned int id, PoseType &x_k_i, double w = 0 );

  /** Destructor */
  ~Particle();

  /**
   *  Set particle pose
   *  \param[in] x_k_i pose
   */
  void setPose( PoseType &x_k_i );

  /**
   *  Get particle pose
   *  \param[out] x_k_i pose
   */
  void getPose( PoseType &x_k_i );

  /**
   *  Get particle pose
   *  \return pointer to pose
   */
  const PoseType* getPose();

  /**
   *  Set particle importance weight
   *  \param[in] w weight
   */
  void setWeight( double w );

  /** 
   *  Get particle importance weight
   *  \return particle weight
   */
  double getWeight();

  /**
   * Get the particle id
   * \return particle id
   */
  unsigned int getId();

  /**
   * Get the particle id of the particle which spawned 
   * the current one after resampling
   * \return parent particle id
   */
  unsigned int getParentId();

  /** 
   * Set the particle id
   * \param[in] id particle id
   */
  void setId( unsigned int id );

  /**
   * Set the id of this particle's parent from resampling
   * \param[in] id parent id
   */
  void setParentId( unsigned int id );

  /** 
   * Copy the state from this particle to another particle
   * \param[out] p particle to which data is copied to
   */
  void copyStateTo( Particle<PoseType>* p);

protected:

  unsigned int id_; /**< particle id number */
  unsigned int idParent_; /** id of particle which this one is spawned from */
  PoseType x_k_i_; /**< robot pose at timestep k **/
  double w_; /**< Particle weight */

};

// Implementation

template< class PoseType >
Particle<PoseType>::Particle(){
  id_ = 0;
  idParent_ = id_;
  x_k_i_ = PoseType();
  w_ = 0;
}

template< class PoseType >
Particle<PoseType>::Particle( unsigned int id, PoseType &x_k_i, double w ){
  id_ = id;
  idParent_ = id_;
  x_k_i_ = x_k_i;
  w_ = w;
}

template< class PoseType >
Particle<PoseType>::~Particle(){}

template< class PoseType >
void Particle<PoseType>::setPose( PoseType &x_k_i ){ 
  x_k_i_ = x_k_i;
}

template< class PoseType >
void Particle<PoseType>::getPose( PoseType &x_k_i ){
  x_k_i = x_k_i_;
}

template< class PoseType >
const PoseType* Particle<PoseType>::getPose(){
  return &x_k_i_;
}

template< class PoseType >
void Particle<PoseType>::setWeight( double w ){ 
  w_ = w;
}

template< class PoseType >
double Particle<PoseType>::getWeight(){ 
  return w_; 
}

template< class PoseType >
unsigned int Particle<PoseType>::getId(){ 
  return id_; 
}

template< class PoseType >
unsigned int Particle<PoseType>::getParentId(){ 
  return idParent_; 
} 

template< class PoseType >
void Particle<PoseType>::setId( unsigned int id ){
  id_ = id;
}

template< class PoseType >
void Particle<PoseType>::setParentId( unsigned int id ){
  idParent_ = id;
}

template< class PoseType >
void Particle<PoseType>::copyStateTo(Particle<PoseType>* p){
  p->x_k_i_ = x_k_i_;
}

#endif
