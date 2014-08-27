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


/** This example shows how lexicographic ordering for pairing landmarks and measurements
 *  can be performed when we have clutter and mis-detections. 
 */

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <utility>
#include "CostMatrix.hpp"
#include "PermutationLexicographic.hpp"

#include <boost/timer/timer.hpp>
#include <boost/format.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

using namespace rfs;

int main(int argc, char *argv[])
{

  printf("5 landmarks and 3 measurements\n");
  printf("First 5 columns indicate the measurement number to which the landmark is associated\n");
  printf("A value == 3 in the first 5 columns imply that a landmark was mis-detected\n");
  printf("A value == 3 in the last 3 columns represent that the measurement is an outlier\n\n");

  unsigned int const nM = 3;
  unsigned int const nZ = 5;
  uint* ordering = new uint[nM + nZ];
 
  PermutationLexicographic pl(nM, nZ, true);
  unsigned int nPerm = pl.next(ordering);
  while( nPerm != 0){
    printf("[%d] %d %d %d | %d %d %d %d %d ", nPerm, ordering[0], ordering[1], ordering[2], ordering[3], ordering[4], ordering[5], ordering[6], ordering[7]);

    printf("| Outliers: ");
    for(int i = nM; i < nM + nZ; i++){
      if(ordering[i] < nZ)
	printf("%d ",ordering[i]);
    }
    printf("\n");
    nPerm = pl.next(ordering);
  }
  
  delete[] ordering;
  return 0;
}
