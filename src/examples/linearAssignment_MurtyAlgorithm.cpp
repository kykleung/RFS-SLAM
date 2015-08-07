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


/** This example shows how Murty's algorithm is used successively select the next-best 
 *  linear assignment when there are clutter and mis-detections
 */

#include <stdio.h>
//#include <vector>
#include <utility>
#include "MurtyAlgorithm.hpp"
#include "BruteForceAssignment.hpp"

#include <boost/timer/timer.hpp>
#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

using namespace rfs;

int main(int argc, char *argv[])
{
  printf("4 landmarks and 3 measurements\n");
  printf("Some measurement log-likelihood matrix:\n");
  printf("Rows 0 - 2 are real landmarks\n");
  printf("Rows 3 - 5 are clutter-generating landmarks\n");
  printf("Cols 0 - 2 are real measurements\n");
  printf("Cols 3 - 5 are miss measurements\n");

  boost::mt19937 rng;
  boost::uniform_01<> ud;
  boost::variate_generator<boost::mt19937, boost::uniform_01<> > gen(rng, ud);

  int const nR1 = 3;
  int const nC1 = 4;
  int const C_size = nR1 + nC1;
  double bigNegNum = -1000;

  double** C = new double*[C_size];
  for(int i = 0; i < C_size; i++){
    C[i] = new double[C_size];
    for(int j = 0; j < C_size; j++){
      if(i < nR1 && j < nC1){
	C[i][j] = log(gen());
      }else if(i < nR1 && j >= nC1){
	if(j - nC1 == i)
	  C[i][j] = log(gen());
	else
	  C[i][j] = bigNegNum;
      }else if(i >= nR1 && j < nC1){
	if(i - nR1 == j)
	  C[i][j] = log(gen());
	else
	  C[i][j] = bigNegNum;
      }else{
	C[i][j] = 0;
      }
    }
  }

  for(int i = 0; i < C_size; i++){
    for(int j = 0; j < C_size; j++){
      printf("%f  ", C[i][j]);
    }
    printf("\n");
  }

  boost::timer::cpu_timer timer;

  Murty murtyAlgo(C, C_size);
  Murty::Assignment a;
  int rank;
  double score;
  murtyAlgo.setRealAssignmentBlock(nR1, nC1);
  
  for(int k = 0; k < 100; k++){
    rank = murtyAlgo.findNextBest(a, score);
    if(rank == -1 || score < bigNegNum)
      break;
    printf("[%d : %f] %d %d %d %d %d %d %d\n", rank, score, a[0], a[1], a[2], a[3], a[4], a[5], a[6]); 
  }  
  boost::timer::cpu_times t = timer.elapsed();
  std::cout << boost::format("Elapsed time: %1% [ns]\n\n") % t.wall;


  // Brute force approach
  printf("Validate result with a brute force approach\n");
  BruteForceLinearAssignment bf;
  int d = 0;
  unsigned int** bfa;
  double* bfs;
  int nbfa = bf.run(C, C_size, bfa, bfs);
  for(int i = 0; i < 300; i++){
    if(bfs[i] < bigNegNum)
      break;
    if(i == 0 || bfs[i] != bfs[i-1]){
      d++;
      printf("[%d : %f] %d %d %d %d %d %d %d\n", d, bfs[i], bfa[i][0], bfa[i][1], bfa[i][2], bfa[i][3], bfa[i][4], bfa[i][5], bfa[i][6]);
    }
  }

  for(int i = 0; i < C_size; i++){
    delete[] C[i];
  }
  delete[] C;

  return 0;
}
