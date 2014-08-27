/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung
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

#ifndef HUNGARIAN_METHOD_HPP
#define HUNGARIAN_METHOD_HPP

#include "CostMatrix.hpp"
#include <queue>
#include <vector>

namespace rfs{

/** 
 * \class HungarianMethod
 * This class is an implementation of the Hungarian method for a linear assignment
 * problem specified by a NxN cost matrix. It has complexity O(N^3).
 * This code is referenced from the Top Coder tutorial on the Hungarian method: 
 * http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=hungarianAlgorithm 
 * The following pdf is also a good reference for the method:
 * www.cse.ust.hk/~golin/COMP572/Notes/Matching.pdf
 * \brief The Hungarian Method for linear assignmnet
 * \author Keith Leung
 */
class HungarianMethod
{
public:

  /** Default constructor */
  HungarianMethod();

  /** Default destructor */
  ~HungarianMethod();

  /**
   * Run the Hungarian method
   * \param[in] C square score / cost matrix.
   * \param[in] n size of cost matrix
   * \param[out] soln assignment solution, memory needs to be allocated by caller 
   * \param[out] cost assignment solution cost, memory needs to be allocated by caller 
   * \param[in] maximize true if we want to find maximum score, false for minimum score
   * \param[in] debug creates debug printouts if true
   * \return whether a solution has been found
   */
  bool run(double** C, int n, int* soln, double* cost, bool maximize = true, bool debug = false );


  /**
   * Run the Hungarian method
   * \param[in] C cost matrix
   * \param[out] soln assignment solution, memory needs to be allocated by caller 
   * \param[out] cost assignment solution cost, memory needs to be allocated by caller 
   * \param[in] maximize true if we want to find maximum score, false for minimum score
   * \return whether a solution has been found
   */
  bool run(CostMatrix &C, int* soln, double* cost, bool maximize = true);

private:


};

}

#endif
