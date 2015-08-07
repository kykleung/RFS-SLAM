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

#include <assert.h>
#include <boost/multi_array.hpp>
#include <cfloat>
#include <cmath>
#include "CostMatrix.hpp"
#include <queue>
#include <stdio.h>
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

    /** Cost Matrix using double pointers */
    typedef double** CostMat_DP;

    /** Cost matrix using boost multi_array*/
    typedef boost::multi_array<double, 2> CostMat;

    /** Default constructor */
    HungarianMethod(){}

    /** Default destructor */
    ~HungarianMethod(){}

    /**
     * Run the Hungarian method
     * \tparam TCostMat Cost matrix type [CostMat_DP or CostMat]
     * \param[in] C square score / cost matrix.
     * \param[in] n size of cost matrix
     * \param[out] soln assignment solution, memory needs to be allocated by caller 
     * \param[out] cost assignment solution cost, memory needs to be allocated by caller 
     * \param[in] maximize true if we want to find maximum score, false for minimum score
     * \param[in] debug creates debug printouts if true
     * \return whether a solution has been found
     */
    template<typename TCostMat>
    bool run(TCostMat &C, int n, int* soln, double* cost, bool maximize = true, bool debug = false );

  };

  /// Hungarian Method Class Implementation ///

  template<typename TCostMat>
  bool HungarianMethod::run(TCostMat &C, int n, int* soln, double* cost, bool maximize, bool debug){

    // Initialize
    double* lx = new double[n]; // label for X
    double* ly = new double[n]; // label for Y
    int* xy = new int[n]; // match to x = X[n] 
    int* yx = new int[n]; // match to y = Y[n]
    bool* S = new bool[n];
    bool* T = new bool[n];
    bool* NS = new bool[n]; // neighbourhood of S, a subset or equal to Y
    double* slack = new double[n];
  
    std::queue<int> q; // queue for breadth-first-search (bfs)
    bool* x_q = new bool[n]; // flag indicating x has already been queued in bfs
    bool* y_q = new bool[n]; // flag indicating y has already been queued in bfs
    int* p = new int[2*n]; // parent of vertices in bfs
  
    int x, x_t, y; // indices
    int root; // root index
    bool pickFreeVertex = true; // flag to pick a new root, step 2
    bool updateLabel; // flag to perform update label, step 3

    for(int x = 0; x < n; x++){
      xy[x] = -1;
      S[x] = false;
    }
    for(int y = 0; y < n; y++){
      yx[y] = -1;
      T[y] = false;
    }

    double offset = 0;
    if(!maximize){
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
    }else{
      for(int x = 0; x < n; x++){
	for(int y = 0; y < n; y++){
	  if(C[x][y] < offset){
	    offset = C[x][y];
	  }
	}
      }
      for(int x = 0; x < n; x++){
	for(int y = 0; y < n; y++){
	  C[x][y] -= offset;
	}
      }
    }
  
    if(debug){
      printf("Score / Cost Matrix:\n");
      for(int x = 0; x < n; x++){
	for(int y = 0; y < n; y++){
	  printf("%f   ", C[x][y]);
	}
	printf("\n");
      }
    }
  

    // Step 1 - initial labeling and matching M in equality graph El
    for( int x = 0; x < n; x++){

      lx[x] = 0.0;
      ly[x] = 0.0;
      for( int y = 0; y < n; y++){
	if( C[x][y] >= lx[x] ){
	  lx[x] = C[x][y];
	  xy[x] = y; // temporary matching of x to y, but there may be conflicts
	}
      }

      // check if y = xy[x] is matched to another x_t = yx[y] = yx[xy[x]]
      int y = xy[x];
      x_t = yx[y];
      //printf("Trying to match x[%d] to y[%d]\n", x, xy[x]); 
      if(yx[y] != -1){ // matched with another x, conflict! Keep the one with the best score
	//printf("  y[%d] is already matched with x[%d]\n", y, yx[y]);
	if( C[x][y] > C[x_t][y] ){
	  xy[x_t] = -1; // current match is better than other match, remove other match
	  yx[y] = x;
	  //printf("  x[%d] is no longer matched with y[%d]\n", x_t, xy[x]);
	  //printf("  x[%d] is now matched with y[%d]\n", x, xy[x]);
	}else{
	  xy[x] = -1; // current match is worse than other match, keep the other match
	  //printf("  x[%d] cannot be matched\n", x);
	}
      }else{
	yx[y] = x;
	//printf("  x[%d] is now matched with y[%d]\n", x, xy[x]);
      }

    }

    if(debug){
      printf("Initial labeling:\n");
      for( int x = 0, y = 0; x < n; x++, y++){
	printf("  lx[%d] = %f     ly[%d] = %f\n", x, lx[x], y, ly[y]);
      }

      printf("Initial matching:\n");
      for( int x = 0; x < n; x++ )
	if( xy[x] != -1)
	  printf("  x[%d] ----- y[%d]\n", x, xy[x]);
	else
	  printf("  x[%d] -----      \n", x);
      for( int y = 0; y < n; y++ )
	if( yx[y] != -1)
	  printf("  y[%d] ----- x[%d]\n", y, yx[y]);
	else
	  printf("  y[%d] -----      \n", y);

      printf("Initial Slack Table:\n");
      for(int x = 0; x < n; x++){
	for(int y = 0; y < n; y++){
	  printf("%e   ", lx[x] + ly[y] - C[x][y]);
	}
	printf("\n");
      }
  
    }

    while(true){

      // Step 2 - Check if we have perfect matching, if not, pick a new root
    
      if(pickFreeVertex){

	for(x = 0; x < n; x++){
	  S[x] = false;
	}
	for(y = 0; y < n; y++){
	  T[y] = false;
	  NS[y] = false;
	}
	for(x = 0; x < n; x++){
	  if(xy[x] == -1)
	    break; // x is not matched and therefore a free vertex
	}
	if(x == n){ // no more free / unmatched x, we have found a solution

	  // undo offset if necessary
	  if(offset != 0){
	    for(x = 0; x < n; x++){
	      for(y = 0; y < n; y++){
		C[x][y] = C[x][y] + offset;
	      }
	    }
	  }

	  // undo negative if necessary
	  if(!maximize){
	    for(x = 0; x < n; x++){
	      for(y = 0; y < n; y++){
		C[x][y] *= -1;
	      }
	    }
	  }

	  *cost = 0;
	  for(x = 0; x < n; x++){
	    soln[x] = xy[x];
	    *cost += C[x][ xy[x] ];
	  }
	  /*printf("Solution found with score = %f\n", *cost);
	    for( int x = 0; x < n; x++ )
	    printf("  x[%d] ----- y[%d]\n", x, xy[x]);*/

	  // Cleanup
	  delete[] lx;
	  delete[] ly;
	  delete[] xy;
	  delete[] yx;
	  delete[] S;
	  delete[] T;
	  delete[] NS;
	  delete[] slack; 
	  delete[] x_q;
	  delete[] y_q;
	  delete[] p;

	  return true;
	}
	root = x; // root of alternating tree
	S[x] = true; // put the free vertex in S
	if(debug)
	  printf("Slack:\n");
	for(y = 0; y < n; y++ ){
	  slack[y] = lx[x] + ly[y] - C[x][y];
	  if(debug)
	    printf("    s[%d] = %f\n", y, slack[y]);
	  if(  fabs(slack[y]) < 1e-14 ){
	    slack[y] = 0;
	    NS[y] = true;
	  }
	}

	if(debug){
	  printf("New root selected\n");  
	  printf("S = { ");
	  for( int x = 0; x < n; x++ ){
	    if(S[x])
	      printf("x[%d] ", x);
	  }
	  printf(" }\n");
	  printf("T = { ");
	  for( int y = 0; y < n; y++ ){
	    if(T[y])
	      printf("y[%d] ", y);
	  }
	  printf(" }\n");
	  printf("NS = { ");
	  for( int y = 0; y < n; y++ ){
	    if(NS[y])
	      printf("y[%d] ", y);
	  }
	  printf(" }\n");
	}
      }

      // Step 3 - Update labels (if neighbourhood N(S) == T)
      updateLabel = true;
      for(y = 0; y < n; y++ ){
	if( NS[y] != T[y] ){
	  updateLabel = false;
	  break;
	}
      }
      if(updateLabel){ // After updating label, N(S) != T


	if(debug){
	  printf("Slack Table:\n");
	  for(int x = 0; x < n; x++){
	    for(int y = 0; y < n; y++){
	      printf("%e   ", lx[x] + ly[y] - C[x][y]);
	    }
	    printf("\n");
	  }
	}

	double a = DBL_MAX;
	//calculate delta using slack    
	for (y = 0; y < n; y++){   	     
	  if (!T[y]){
	    a = fmin(a, slack[y]);
	  }
	}

	if(debug){
	  printf("Update labels:\n");
	  printf("  minimum slack of !T = %f\n",a);
	}
    
	//update X labels
	for (x = 0; x < n; x++){            
	  if (S[x]){ 
	    lx[x] -= a;
	  }
	}
    
	//update Y labels
	for (y = 0; y < n; y++){            
	  if (T[y]){ 
	    ly[y] += a;
	  }
	}

	//update slack
	for (y = 0; y < n; y++){            
	  if (!T[y]){
	    slack[y] -= a;
	  }
	  if(slack[y] == 0){
	    NS[y] = true;
	  }
	}

	if(debug){
	  for( int x = 0, y = 0; x < n; x++, y++)
	    printf("  lx[%d] = %f     ly[%d] = %f\n", x, lx[x], y, ly[y]);
	  printf("  Update slack variables:\n");
	  for( int y = 0; y < n; y++ )
	    printf("    s[%d] = %f\n", y, slack[y]);

	  printf("Slack Table:\n");
	  for(int x = 0; x < n; x++){
	    for(int y = 0; y < n; y++){
	      printf("%e   ", lx[x] + ly[y] - C[x][y]);
	    }
	    printf("\n");
	  }
	}

      }//end updateLabel

      // Step 4 - N(S) != T, select y in N(S) but not in T
      //printf("Finding y in NS and not T\n");
      for(y = 0; y < n; y++ ){
	if( NS[y] && !T[y] )
	  break;
      }
      x_t = yx[y]; // match to y
      if(x_t == -1){ // if y is free, path root -> y is an augmenting path, goto step 2

	bool augmentingPathFound = false;
      
	if(debug){
	  printf("  Found free y[%d]\n", y);
	}

	// Breadth first search for path root -> y
	// We need to make distinct indices for vertices in X and in Y
	// So in this search only, indices for X are 0, 1, 2, 3 ... n-1
	// and indices for Y are n+0, n+1, n+2 ... 2n-1 when storing indices in the queue
	int target = y + n;
	std::queue<int> q_empty;
	std::swap(q, q_empty); // clear the queue
	q.push(root);
	for(x = 0; x < n; x++){
	  x_q[x] = false;
	  y_q[x] = false;
	}
	x_q[root] = true;
	for(x = 0; x < 2*n; x++){
	  p[x] = -1;
	}
	int t = q.front();
	while(!q.empty()){
 
	  t = q.front();
	  if(debug){
	    if( t < n){
	      printf("Breadth first search queue front: x[%d]\n",t);
	    }else{
	      printf("Breadth first search queue front: y[%d]\n",t-n);	  
	    }
	  }

	  if (t == target){
	    // path found, perform augmentation in inverse the matching along this path
	    if(debug){
	      printf("  Augmentation path found:\n    y[%d]", t-n);
	    }
	    while(t != root){
	      if( t >= n){ // t is in Y, p[t]--t should be a match
		x_t = p[t];
		xy[x_t] = t-n;
		yx[t-n] = x_t;
	      }
	      t = p[t];

	      if(debug){
		if( t < n ){
		  printf(" --- x[%d]", t);
		}else{
		  printf(" --- y[%d]", t-n);
		}
	      }
	    }

	    if(debug){
	      printf("\nCurrent pairing\n");
	  
	      for( int x = 0; x < n; x++ )
		if( xy[x] != -1){
		  printf("  x[%d] ----- y[%d]\n", x, xy[x]);
		}else{
		  printf("  x[%d] -----      \n", x);
		}
	      for( int y = 0; y < n; y++ )
		if( yx[y] != -1){
		  printf("  y[%d] ----- x[%d]\n", y, yx[y]);
		}else{
		  printf("  y[%d] -----      \n", y);
		}
	    }
	    augmentingPathFound = true;
	    break;
	  }
	  q.pop();

	  if( t < n){ // t is a vertex in X 
	    for(y = 0; y < n; y++){ // look through children of t
	      // check if y is in equality graph and not already queued and unmatched
	      if(fabs(lx[t] + ly[y] - C[t][y]) < 1e-12 && !y_q[y] && xy[t] != y){
		if(debug){
		  printf("Breadth first search inserting y[%d] in queue with parent x[%d]\n", y, t); 
		}
		y_q[y] = true;
		p[y+n] = t;
		q.push( y+n );
	      }
	    }
	  }else{ // t >= n and is a vertex in Y
	    t -= n; // fix the index back to original
	    for(x = 0; x < n; x++){ // look through children of t
	      // check if x is in equality graph and in set S and not already queued and matched
	      if(fabs(lx[x] + ly[t] - C[x][t]) < 1e-12 && S[x] && !x_q[x] && yx[t] == x){
		if(debug){
		  printf("Breadth first search inserting x[%d] in queue with parent y[%d]\n", x, t);
		}
		x_q[x] = true;
		p[x] = t+n;
		q.push( x );
	      }
	    }
	  }

	} // end while

	if(!augmentingPathFound){
	  printf("Cannot find alternating path\n");
	  /*printf("Slack table:\n");
	    for(int x = 0; x < n; x++){
	    for(int y = 0; y < n; y++){
	    printf("%e   ", lx[x] + ly[y] - C[x][y]);
	    }
	    printf("\n");
	    }*/
	  return false;
	}

	pickFreeVertex = true; // go back to step 2

      }else{ // y is matched to x_t, add to alternating tree, go to step 3

	if(debug){
	  printf("  Found y[%d] matched with x[%d]\n", y, x_t);
	  printf("  Examine slack for x[%d]:\n", x_t);
	}
	S[x_t] = true;
	T[y] = true;
	for(int y_t = 0; y_t < n; y_t++ ){
	  if(debug)
	    printf("    s[%d -- %d] = %f\n", x_t, y_t, lx[x_t] + ly[y_t] - C[x_t][y_t]);
	  if( fabs(lx[x_t] + ly[y_t] - C[x_t][y_t]) < 1e-14 )
	    NS[y_t] = true;
	}
      
	if(debug){
	  printf("  S = { ");
	  for( int x = 0; x < n; x++ ){
	    if(S[x]){
	      printf("x[%d] ", x);
	    }
	  }
	  printf(" }\n");
	  printf("  T = { ");
	  for( int y = 0; y < n; y++ ){
	    if(T[y]){
	      printf("y[%d] ", y);
	    }
	  }
	  printf(" }\n");
	  printf("  NS = { ");
	  for( int y = 0; y < n; y++ ){
	    if(NS[y]){
	      printf("y[%d] ", y);
	    }
	  }
	  printf(" }\n");
	}

	// It is necessary to update the slack variables here since S now includes x_t
	for(int y = 0; y < n; y++ ){
	  double slack_x_t = lx[x_t] + ly[y] - C[x_t][y]; 
	  if( slack_x_t < slack[y]){
	    slack[y] = slack_x_t;
	  }
	}

	if(debug){
	  printf("  Update slack variables for !T:\n");
	  for( int y = 0; y < n; y++ )
	    if(!T[y])
	      printf("    s[%d] = %f\n", y, slack[y]);
	}      

	pickFreeVertex = false; // go back to step 3

      }

    }

  }

  template<>
  bool HungarianMethod::run(CostMatrix &C, int n, int* soln, double* cost, bool maximize, bool debug);

}

#endif
