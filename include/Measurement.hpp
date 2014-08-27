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

// Measurement Classes

#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include "RandomVec.hpp"

namespace rfs{

/** Definition for a NULL measurement */
typedef RandomVec < ::Eigen::Matrix<double, 1, 1>,
		    ::Eigen::Matrix<double, 1, 1> > NullInput;

/** Definition for 1d measurement */
typedef RandomVec < ::Eigen::Matrix<double, 1, 1>,
		    ::Eigen::Matrix<double, 1, 1> > Measurement1d;

/** Definition for 2d measurement */
typedef RandomVec < ::Eigen::Vector2d, ::Eigen::Matrix2d> Measurement2d;

/** Definition for 3d measurement */
typedef RandomVec < ::Eigen::Vector3d, ::Eigen::Matrix3d> Measurement3d;


/** Definition for 1d odometry */
typedef RandomVec <  ::Eigen::Matrix<double, 1, 1>,
		     ::Eigen::Matrix<double, 1, 1> > Odometry1d;


/********** 2d odometry measurement **********/

/**
 * \class Odometry2d
 * \brief A class for 2d odometry measurements for a 2d motion model
 * \author Keith Leung
 * \note This class is derived from so that we can make a customized constructor for convenience.
 */
class Odometry2d : public RandomVec< ::Eigen::Vector3d, ::Eigen::Matrix3d >
{
public:

  /** For using Eigen fixed-size matrices with STL containers */
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  
  /** Default constructor */
  Odometry2d();
 
  /** 
   * Constructor - defined only for our convenience
   * \param[in] x measurement vector
   * \param[in] Sx measurement uncertainty covariance
   * \param[in] t time
   */
  Odometry2d(Vec &x, Mat &Sx, TimeStamp t = TimeStamp());

  /** 
   * Constructor - defined only for our convenience
   * \param[in] x measurement vector
   * \param[in] t time
   */
  Odometry2d(Vec &x, TimeStamp t = TimeStamp());


  /** 
   * Constructor - defined only for our convenience
   * \param[in] dx_k_km x-displacement from frame k-1
   * \param[in] dy_k_km y-displacement from frame k-1
   * \param[in] dtheta_k_km rotational displacement from frame k-1
   * \param[in] vardx_k_km variance in dx_k_km
   * \param[in] vardy_k_km variance in dy_k_km
   * \param[in] vartheta_k_km variance in dtheta_k_km
   * \param[in] t time of odometry reading
   */
  Odometry2d(double dx_k_km, double dy_k_km, double dtheta_k_km,
	     double vardx_k_km, double vardy_k_km, double vartheta_k_km,
	     TimeStamp t = TimeStamp());

  /** Destructor */
  ~Odometry2d();
};

}

#endif
