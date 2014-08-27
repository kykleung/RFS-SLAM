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

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <cstddef>
#include "RandomVec.hpp"

namespace rfs{

/** 
 *  \class Particle
 *  \brief A class for a particle for the particle filter
 *  \tparam PoseType RandomVec derived class to represent the state of a particle
 *  \tparam DataType Some user-defined class for carrying extra data / information with a particle
 *  \author Keith Leung
 */
template< class PoseType, class DataType = int>
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
  void getPose( PoseType &x_k_i ) const;

  /**
   *  Get particle pose
   *  \return pointer to pose
   */
  const PoseType* getPose() const;

  /**
   *  Set particle importance weight
   *  \param[in] w weight
   */
  void setWeight( double w );

  /** 
   *  Get particle importance weight
   *  \return particle weight
   */
  double getWeight() const;

  /**
   * Get the particle id
   * \return particle id
   */
  unsigned int getId() const;

  /**
   * Get the particle id of the particle which spawned 
   * the current one after resampling
   * \return parent particle id
   */
  unsigned int getParentId() const ;

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
  void copyStateTo( Particle<PoseType, DataType>* p);

  /** 
   * Copy the extra data from this particle to another particle
   * \warning DataType may need a non-default copy constructor
   * \param[out] p particle to which data is copied to
   */
  void copyDataTo( Particle<PoseType, DataType>* p);

  /** 
   * Get the extra data pointer 
   * \return pointer to the data
   */
  DataType* getData();

  /**
   * Set the extra data pointer
   * \param[in] dataPtr pointer to the data
   */
  void setData(DataType* dataPtr);

  /**
   * Delete the extra data
   */
  void deleteData();

protected:

  unsigned int id_; /**< particle id number */
  unsigned int idParent_; /** id of particle which this one is spawned from */
  PoseType x_k_i_; /**< robot pose at timestep k **/
  double w_; /**< Particle weight */
  DataType* data_; /**< Pointer to extra data */

};

// Implementation

template< class PoseType, class DataType >
Particle<PoseType, DataType>::Particle(){
  id_ = 0;
  idParent_ = id_;
  x_k_i_ = PoseType();
  w_ = 0;
  data_ = NULL;
}

template< class PoseType, class DataType >
Particle<PoseType, DataType>::Particle( unsigned int id, PoseType &x_k_i, double w ){
  id_ = id;
  idParent_ = id_;
  x_k_i_ = x_k_i;
  w_ = w;
  data_ = NULL;
}

template< class PoseType, class DataType >
Particle<PoseType, DataType>::~Particle(){
  deleteData();
}

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::setPose( PoseType &x_k_i ){ 
  x_k_i_ = x_k_i;
}

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::getPose( PoseType &x_k_i ) const{
  x_k_i = x_k_i_;
}

template< class PoseType, class DataType >
const PoseType* Particle<PoseType, DataType>::getPose() const{
  return &x_k_i_;
}

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::setWeight( double w ){ 
  w_ = w;
}

template< class PoseType, class DataType >
double Particle<PoseType, DataType>::getWeight() const{ 
  return w_; 
}

template< class PoseType, class DataType >
unsigned int Particle<PoseType, DataType>::getId() const{ 
  return id_; 
}

template< class PoseType, class DataType >
unsigned int Particle<PoseType, DataType>::getParentId() const{ 
  return idParent_; 
} 

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::setId( unsigned int id ){
  id_ = id;
}

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::setParentId( unsigned int id ){
  idParent_ = id;
}

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::copyStateTo(Particle<PoseType, DataType>* p){
  p->x_k_i_ = x_k_i_;
}

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::copyDataTo(Particle<PoseType, DataType>* p){
  if( data_ != NULL )
    p->data_ = new DataType( *data_ );
}


template< class PoseType, class DataType >
DataType* Particle<PoseType, DataType>::getData(){
  return data_;
}

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::setData(DataType* dataPtr){
  data_ = dataPtr;
}

template< class PoseType, class DataType >
void Particle<PoseType, DataType>::deleteData(){
  if( data_ != NULL ){
    delete data_;
    data_ = NULL;
  }
}
  
}  

#endif
