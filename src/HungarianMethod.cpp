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

#include "HungarianMethod.hpp"

namespace rfs
{

  template<>
  bool HungarianMethod::run<CostMatrix>(CostMatrix &C, int n, int* soln, double* cost, bool maximize, bool debug){

    double** C_full;
    double** C_reduced;
    unsigned int nRowsFull, nColsFull;
    C.getCostMatrix(C_full, nRowsFull, nColsFull);
    *cost = 0;
    double score_fixed = 0;
    int* i_remap = NULL;
    int* j_remap = NULL;
    int n_reduced = C.getCostMatrixReduced(C_reduced, soln, &score_fixed, i_remap, j_remap);
    if(n_reduced > 1){
      // reduced cost matrix available
      int soln_reduced[n_reduced];
      double score_reduced = 0;
      if( run(C_reduced, n_reduced, soln_reduced, &score_reduced, maximize) ){
	for(int n = 0; n < n_reduced; n++ ){
	  soln[ i_remap[n] ] = j_remap [ soln_reduced[n] ];
	}
	*cost = score_fixed + score_reduced;
      }else{
	return false;
      }
    }else if(n_reduced == 0){
      // solution is already available through reduction
      *cost = score_fixed;
      return true;
    }else{
      // no reduction on C was performed
      return run(C_full, nRowsFull, soln, cost, maximize);
    }
    return false;
  }



}

