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

#ifndef PROCESSMODEL_ODOMETRY1D_HPP
#define PROCESSMODEL_ODOMETRY1D_HPP

#include "ProcessModel.hpp"
#include "Measurement.hpp"
#include "Pose.hpp"

namespace rfs{

/**
 * \class MotionModel_Odometry1d
 * A 1d odometry motion model with translational displacement
 * The 1d model is as follows:
 * \f[ x_{k} = x_{k-1} + \delta x_k\f]
 * \brief A 1d odometry motion model with translational displacement
 * \note Currently the updated state from step does not contain valid
 * covariance information because it is not needed by the RBPHDFilter
 * \author Keith Leung
 */
class MotionModel_Odometry1d : public ProcessModel< Pose1d, Odometry1d >
{
public:

  /** Default constructor */
  MotionModel_Odometry1d(){}

  /** Constructor with process noise input 
   * \param[in] Q additive zero-mean white Gaussian noise variance
   */
  MotionModel_Odometry1d( Pose1d::Mat &Q );

  /** Default destructor */
  ~MotionModel_Odometry1d(){}

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
	     TimeStamp const &dT);

};

}

#endif
