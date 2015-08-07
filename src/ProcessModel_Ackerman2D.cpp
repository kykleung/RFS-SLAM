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

#include "ProcessModel_Ackerman2D.hpp"

using namespace rfs;

MotionModel_Ackerman2d::MotionModel_Ackerman2d(Pose2d::Cov Q): ProcessModel(Q){}

MotionModel_Ackerman2d::MotionModel_Ackerman2d(double h, double l, double dx, double dy, Pose2d::Cov Q): 
h_(h), l_(l), poi_offset_x_(dx), poi_offset_y_(dy), ProcessModel(Q){}

MotionModel_Ackerman2d::~MotionModel_Ackerman2d(){}

void MotionModel_Ackerman2d::setAckermanParams(double h, double l, double dx, double dy){
  h_ = h;
  l_ = l;
  poi_offset_x_ = dx;
  poi_offset_y_ = dy;
}

void MotionModel_Ackerman2d::step( Pose2d &s_k, Pose2d &s_km, AckermanInput &input_k, TimeStamp const &dT){

  double u_v = input_k[0]; // velocity
  double u_r = input_k[1]; // steering
  double dt = dT.getTimeAsDouble();

  Pose2d::Vec x_km = s_km.get();
  double r_km = x_km[2];
  double cos_r_km = cos(r_km);
  double sin_r_km = sin(r_km);
  double tan_u_r = tan(u_r);
  TimeStamp t_km = s_km.getTime();

  u_v = u_v / (1 - tan_u_r * h_ / l_);

  Pose2d::Vec dx;
  dx << dt * ( u_v * cos_r_km - u_v / l_ * tan_u_r * (poi_offset_x_ * sin_r_km + poi_offset_y_ * cos_r_km ) ),
        dt * ( u_v * sin_r_km + u_v / l_ * tan_u_r * (poi_offset_x_ * cos_r_km - poi_offset_y_ * sin_r_km ) ),
        dt * u_v / l_ * tan_u_r;
  Pose2d::Vec x_k = x_km + dx;
  if(x_k(2) > PI){
    x_k(2) -= 2 * PI;
  }else if(x_k(2) < -PI){
    x_k(2) += 2 * PI;
  } 

  s_k.set(x_k, t_km + dT);

}
