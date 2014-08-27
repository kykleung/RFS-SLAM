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

#include "PermutationLexicographic.hpp"
#include <algorithm>
#include <stdio.h>

namespace rfs
{

PermutationLexicographic::PermutationLexicographic(unsigned int nM, 
						   unsigned int nZ, 
						   bool includeClutter) : nM_(nM), nZ_(nZ), nP_(0), last_(false)
{

  if(nM != nZ)
    includeClutter = true;
  
  oSize_ = 0;
  if(includeClutter == false){
    oSize_ = nM;
  }else{
    oSize_ = nM + nZ;
  }

  o_ = new unsigned int[oSize_];
  for(int i = 0; i < oSize_; i++){
    if( i < nZ )
      o_[i] = i;
    else
      o_[i] = nZ;
  }
  
  
}

PermutationLexicographic::~PermutationLexicographic(){

  delete[] o_;

}

unsigned int PermutationLexicographic::next(unsigned int* permutation){

  if(last_){
    for(int i = 0; i < oSize_; i++){
      permutation[i] = 0;
    } 
    return 0;
  }

  for(int i = 0; i < oSize_; i++){
    permutation[i] = o_[i];
  }
    
  unsigned int u = nM_;
  unsigned int v = oSize_ - 1;

  while(u < v){ // Fast forward permutations to the next meaningful one
    std::swap(o_[u], o_[v]);
    u++;
    v--;
  }
  
  last_ = !std::next_permutation(o_, o_ + oSize_);
  nP_++;
  return nP_;

}

}
