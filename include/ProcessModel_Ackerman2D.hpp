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

#ifndef PROCESSMODEL_ACKERMAN2D_HPP
#define PROCESSMODEL_ACKERMAN2D_HPP

#include "ProcessModel.hpp"
#include "Measurement.hpp"
#include "Pose.hpp"

namespace rfs{

  typedef Pose<2, 1, 1> AckermanInput;

  /**
   * \class MotionModel_Ackerman2d
   * \brief Ackerman motion model.
   *
   * The inputs are rear wheels average velocity and front wheel steering angle.
   * The vehicle center is defined at the center between the rear wheels.
   * The point of interest that we wish to move can be offset from the vehicle center.
   * Source: http://www-personal.acfr.usyd.edu.au/nebot/experimental_data/modeling_info/Ute_modeling_info.htm
   * \note Currently the updated state from step does not contain valid
   * covariance information because it is not needed by the RBPHDFilter
   * \author Keith Leung
   */
  class MotionModel_Ackerman2d : public ProcessModel< Pose2d, AckermanInput >
  {
  public:

    /**
     * \brief Constructor.
     * \param[in] Q additive zero-mean white Gaussian noise (added to the output)
     */
    MotionModel_Ackerman2d(Pose2d::Cov Q = Pose2d::Cov::Zero());

    /** 
     * \brief Constructor.
     * \param[in] h Rear wheel lateral offset from vehicle centerline
     * \param[in] l Distance between front and rear wheels
     * \param[in] dx x offset of the point of interest from the vehicle center
     * \param[in] dy y offset of the point of interest from the vehicle center
     * \param[in] Q additive zero-mean white Gaussian noise (added to the output)
     */
    MotionModel_Ackerman2d(double h, double l, double dx, double dy, Pose2d::Cov Q = Pose2d::Cov::Zero());

    /** Default destructor */
    ~MotionModel_Ackerman2d();

    /**
     * \brief Set the parameters of the Ackerman model
     * \param[in] h Rear wheel lateral offset from vehicle centerline
     * \param[in] l Distance between front and rear wheels
     * \param[in] dx x offset of the point of interest from the vehicle center
     * \param[in] dy y offset of the point of interest from the vehicle center
     */
    void setAckermanParams(double h, double l, double dx, double dy);

    /** 
     * This overrides the virtual function in the parent class for
     * determining the pose at time-step k from pose at time-step k-1.
     * \note Currently the updated state from step does not contain valid
     * covariance information because it is not needed by the RBPHDFilter
     * \param[out] s_k pose at current time-step k
     * \param[in] s_km pose at previous time-step k-1
     * \param[in] input_k input to process model
     * \param[in] dT size of time-step
     */
    void step( Pose2d &s_k, Pose2d &s_km, AckermanInput &input_k, TimeStamp const &dT);

  private:
    
    double h_; /**< Rear wheel lateral offset from vehicle centerline */
    double l_; /**< Distance between front and rear wheels */
    double poi_offset_x_; /**< x offset of the point of interest from the vehicle center */
    double poi_offset_y_; /**< y offset of the point of interest from the vehicle center */
  
  };

}

#endif
