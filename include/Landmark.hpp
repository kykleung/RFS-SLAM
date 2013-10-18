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

// Landmark classes for defining map feature state
// Keith Leung 2013

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include "RandomVec.hpp"

/** 
 * \class Landmark
 * \brief An abstract class for defining landmark state
 * \author Keith Leung
 */

template<class VecType, class MatType, class DescriptorType = int>
class Landmark : public RandomVec<VecType, MatType>
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** Default constructor */
  Landmark(){}

  /** 
   * Constructor
   */
  Landmark(VecType &x, MatType &Sx){
    this->set(x, Sx);
  }

  /** 
   * Constructor
   */
  Landmark(VecType &x, MatType &Sx, DescriptorType &d){
    this->set(x, Sx);
    desc_ = d;
  }

  /** Default destructor */
  ~Landmark(){};

  /** Set descriptor for landmark 
   *  \param[in] d descriptor
   */
  void setDescriptor(DescriptorType &d){
    desc_ = d;
  }

  /** Get descriptor for landmark 
   *  \param[in] d descriptor
   */
  void getDescriptor(DescriptorType &d){
    d = desc_;
  }

private:
  
  DescriptorType desc_;

};

typedef Landmark< Eigen::Matrix<double, 1, 1>, Eigen::Matrix<double, 1, 1> >
Landmark1d;

typedef Landmark<Eigen::Vector2d, Eigen::Matrix2d> Landmark2d;

#endif
