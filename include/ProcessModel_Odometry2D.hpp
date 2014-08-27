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

#ifndef PROCESSMODEL_ODOMETRY2D_HPP
#define PROCESSMODEL_ODOMETRY2D_HPP

#include "ProcessModel.hpp"
#include "Measurement.hpp"
#include "Pose.hpp"

namespace rfs{

/************* 2d Odometry Motion Model *************/

/**
 * \class MotionModel_Odometry2d
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
class MotionModel_Odometry2d : public ProcessModel< Pose2d, Odometry2d >
{
public:

  /** Default constructor */
  MotionModel_Odometry2d();

  /** Constructor with process noise input 
   * \param Q additive zero-mean white Gaussian noise covariance matrix
   */
  MotionModel_Odometry2d( Pose2d::Mat &Q );

  /** Default destructor */
  ~MotionModel_Odometry2d();

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
   * \param[in] dT size of time-step
   */
  void step( Pose2d &s_k, Pose2d &s_km, Odometry2d &input_k, 
	     TimeStamp const &dT);
  
};

}

#endif
