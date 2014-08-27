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
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "BruteForceAssignment.hpp"
#include "PermutationLexicographic.hpp"

#include <cmath>
#include <stdio.h>

namespace rfs
{

BruteForceLinearAssignment::BruteForceLinearAssignment() : a_(NULL), s_(NULL), nAssignments_(0){}

BruteForceLinearAssignment::~BruteForceLinearAssignment(){
  if(a_ != NULL ){
    for( unsigned int m = 0; m < nAssignments_; m++ ){
      delete[] a_[m];
    }
    delete[] a_;
  }
  if(s_ != NULL ){
    delete[] s_;
  }
}

unsigned int BruteForceLinearAssignment::run(double** C, int n, unsigned int** &a, double* &s, bool maxToMin){

  if(a_ != NULL ){
    for( unsigned int m = 0; m < nAssignments_; m++ ){
      delete[] a_[m];
    }
    delete[] a_;
  }
  if(s_ != NULL ){
    delete[] s_;
  }

  double offset = 0;
  if(!maxToMin){
    for(int x = 0; x < n; x++){
      for(int y = 0; y < n; y++){
	C[x][y] *= -1;
	if(C[x][y] < offset)
	  offset = C[x][y];
      }
    }
    for(int x = 0; x < n; x++){
      for(int y = 0; y < n; y++){
	C[x][y] -= offset;
      }
    }
  }

  nAssignments_ = 0;
  PermutationLexicographic lex(n, n, false);
  while(true){
    
    unsigned int* o = new unsigned int[n];
    int np = lex.next(o);
    if(np == 0){
      delete[] o;
      break;
    }
    nAssignments_ = np;
    assignment as;
    as.a = o;
    as.score = 0;
    if( maxToMin ){
      for(int i = 0; i < n; i++ ){
	as.score += C[i][o[i]];
      }
    }else{
      for(int i = 0; i < n; i++ ){
	as.score -= ( C[i][o[i]] + offset );
      }
    }
    
    pq_.push(as);

  }

  a_ = new unsigned int* [nAssignments_];
  s_ = new double [nAssignments_];
  for(int m = 0; m < nAssignments_; m++){
    assignment as = pq_.top();
    pq_.pop();
    a_[m] = as.a;
    s_[m] = as.score;
  }
  a = a_;
  s = s_;
  return nAssignments_;
}

}
