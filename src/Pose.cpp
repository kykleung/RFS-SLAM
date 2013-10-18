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

#include "Pose.hpp"

/********** Implementation of 1d vechile position state ***********/

Pose1d::Pose1d(){}

Pose1d::Pose1d(double x, double Sx, double t){
  Vec state;
  Mat cov;
  state << x;
  cov << Sx;
  set(state, cov, t);
}

Pose1d::Pose1d(Eigen::Matrix<double, 1, 1> &x, Eigen::Matrix<double, 1, 1> &Sx, double t):
  RandomVec< Eigen::Matrix<double, 1, 1>, Eigen::Matrix<double, 1, 1> >(x, Sx, t){}

Pose1d::Pose1d(double x, double t){
  Vec state;
  state << x;
  set(state, t);
}

Pose1d::Pose1d(Eigen::Matrix<double, 1, 1> &x, double t):
  RandomVec< Eigen::Matrix<double, 1, 1>, Eigen::Matrix<double, 1, 1> >(x, t){}


/********** Implementation of 2d vehicle pose state **********/

Pose2d::Pose2d(){}

Pose2d::Pose2d(Vec &x, Mat &Sx, double t) :
  RandomVec< Eigen::Vector3d, Eigen::Matrix3d >(x, Sx, t){}

Pose2d::Pose2d(Vec &x, double t) : 
  RandomVec< Eigen::Vector3d, Eigen::Matrix3d >(x, t){}

Pose2d::Pose2d( double x, double y, double theta, 
		double var_x, double var_y, double var_theta,
		double t ){
  Vec state;
  state << x, y, theta;
  Mat cov;
  cov << 
    var_x, 0, 0,
    0, var_y, 0,
    0, 0, var_x;
  set(state, cov, t);
}

Pose2d::~Pose2d(){}

