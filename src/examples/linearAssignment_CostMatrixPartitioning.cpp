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


/** This example shows how a cost matrix (or a matrix of measurement likelihoods) can
 *    be partitioned into smaller parts to improve computational efficiency
 */

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <utility>
#include "CostMatrix.hpp"

#include <boost/timer/timer.hpp>
#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

using namespace rfs;

int main(int argc, char *argv[])
{

  boost::mt19937 rng;
  boost::uniform_01<> ud;
  boost::variate_generator<boost::mt19937, boost::uniform_01<> > gen(rng, ud);

  int const C_size = 7;
  double** C;
  CostMatrix likelihoodMat(C, C_size);

  for(int i = 0; i < C_size; i++){
    for(int j = 0; j < C_size; j++){
      C[i][j] = 0;
    }
  }
  C[1][4] = gen();
  C[2][1] = gen();
  C[3][3] = gen();
  C[3][5] = gen();
  C[5][2] = gen();
  C[5][3] = gen();
  C[5][5] = gen();

  printf("Cost matrix:\n");
  for(int i = 0; i < C_size; i++){
    for(int j = 0; j < C_size; j++){
      printf("%f  ", C[i][j]);
    }
    printf("\n");
  }


  boost::timer::cpu_timer timer;

  int nP = likelihoodMat.partition();
  
  boost::timer::cpu_times t = timer.elapsed();
  std::cout << boost::format("Elapsed time: %1% [ns]\n") % t.wall;

  for(int n = 0; n < nP; n++){
    
    printf("Partition %d\n", n);

    unsigned int nCols, nRows;
    double** Cp; 
    unsigned int* rowIdx;
    unsigned int* colIdx;
    bool isZeroPartition = !likelihoodMat.getPartition(n, Cp, nRows, nCols, rowIdx, colIdx);
 
    if(isZeroPartition)
      printf("Zero partition\n");
    printf("Original rows: ");
    for(int i = 0; i < nRows; i++){
      printf("%d  ", rowIdx[i]);
    }
    printf("\n");
    printf("Original cols: ");
    for(int j = 0; j < nCols; j++){
      printf("%d  ", colIdx[j]);
    }
    printf("\n");

    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
	printf("%f  ", Cp[i][j]);
      }
      printf("\n");
    }

  }

  return 0;
}
