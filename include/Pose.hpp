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

// Pose classes for defining vehicle state
// Keith Leung 2013

#ifndef POSE_HPP
#define POSE_HPP

#include "RandomVec.hpp"

namespace rfs{


/********** Define a 1d vehicle pose state **********/

/**
 * \class Pose1d
 * Vehicle position in 1d space
 * \brief Vehicle position
 * \author Keith Leung
 * \note The custom constructors are defined for convenience.
 */
class Pose1d : public RandomVec< ::Eigen::Matrix<double, 1, 1>,
				 ::Eigen::Matrix<double, 1, 1> >
{
public:
  
  /** For using Eigen fixed-size matrices with STL containers */
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** Default constructor */
  Pose1d();

  /** Constructor 
   *  \param[in] x position
   *  \param[in] Sx variance
   *  \param[in] t time
   */
  Pose1d( double x, double Sx, const TimeStamp &t);

  /** Constructor 
   *  \param[in] x position
   *  \param[in] Sx variance
   *  \param[in] t time
   */
  Pose1d(const ::Eigen::Matrix<double, 1, 1> &x, const ::Eigen::Matrix<double, 1, 1> &Sx, const TimeStamp &t);

  /** Constructor 
   *  \param[in] x position
   *  \param[in] t time
   */
  Pose1d( double x, const TimeStamp &t);

  /** Constructor 
   *  \param[in] x position
   *  \param[in] t time
   */
  Pose1d(const ::Eigen::Matrix<double, 1, 1> &x, const TimeStamp &t);

  /** Destructor */
  ~Pose1d(){}

};

/********** 2d vehicle pose state **********/

/**
 * \class Pose2d
 * \f[ \mathbf{x} = \begin{bmatrix} x\\ y\\ \theta \end{bmatrix} \f]
 * \brief 2d vehicle Pose for (x,y) coordinate and rotation
 * \author Keith Leung
 * \note This class is derived from pose only so that we can add a 
 *  custom custructor for our convenience.
 */
class Pose2d : public RandomVec< ::Eigen::Vector3d, ::Eigen::Matrix3d >
{

public:

  /** For using Eigen fixed-size matrices with STL containers */
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** 
   * Default constructor
   */
  Pose2d();

  /** 
   * Constructor - defined only for our convenience
   * \param[in] x pose vector
   * \param[in] Sx pose uncertainty covariance
   * \param[in] t time
   */
  Pose2d(const Vec &x, const Mat &Sx, const TimeStamp &t);

  /** 
   * Constructor - defined only for our convenience
   * \param[in] x pose vector 
   * \param[in] t time
   */
  Pose2d(const Vec &x, const TimeStamp &t);

  /**
   * Constructor - defined only for our convenience
   * \param[in] x x-position
   * \param[in] y y-position
   * \param[in] theta orientation
   * \param[in] var_x x-position variance
   * \param[in] var_y y-position variance
   * \param[in] var_theta theta orientation variance
   * \param[in] t time
   */
  Pose2d( double x, double y, double theta,
	  double var_x, double var_y, double var_theta,
	  const TimeStamp &t );

  /** Default destructor */
  ~Pose2d();

};

}

#endif
