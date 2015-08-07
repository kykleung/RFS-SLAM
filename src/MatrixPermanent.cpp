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

#include "MatrixPermanent.hpp"
#include <float.h>
#include <iostream>

namespace rfs{

  MatPerm::MatPerm(){}

  MatPerm::~MatPerm(){}

  double MatPerm::calc(Eigen::MatrixXd &A){

    if(A.rows() != A.cols()){
      return -DBL_MAX;
    }

    int n = A.rows();
    double x[n];
    double p = 0;
    double s = -1;
    int g[n];
    for(int i = 0; i < n; i++){
      g[i] = 0;
    }

    for(int i = 0; i < n; i++){      
      double row_i_sum = 0;
      for(int j = 0; j < n; j++){
	row_i_sum += A(i,j);
      }
      x[i] = A(i,n-1) - 0.5 * row_i_sum;
    }

    p = s;
    for(int i = 0; i < n; i++){
      p *= x[i];
    }

    for(int k = 2; k <= pow(2, n-1); k++){

      int j = 0;

      if( k % 2 == 0){ // k is even, flip first bit

	j = 0;

      }else{ // k is odd, find first bit that is equal to 1, and flip the next bit

	j = 1;
	while( g[j-1] == 0 ){
	  j++;
	}

      }

      s *= -1;
      double z = 1 - 2 * g[j]; // add or subtract flag for A(i,j)
      g[j] = !g[j]; // flip bit     
      
      /*std::cout << "bits: ";
      for(int i = 0; i < n; i++){
	std::cout << g[i] << " ";
      }
      std::cout << "(z = " << z << ")" << std::endl;
      */

      double x_prod = 1;
      for(int i = 0; i < n; i++){
	x[i] += z * A(i,j);
	x_prod *= x[i];
      }

      p += s * x_prod;
      //std::cout << "p = " << p << std::endl;

    }

    double retval = 2 * p;
    if( n % 2 != 0){
      retval *= -1;
    }
    return retval;
  }

}
