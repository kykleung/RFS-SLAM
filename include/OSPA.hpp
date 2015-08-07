/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2014, Keith Leung
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

#ifndef OSPA_HPP
#define OSPA_HPP

#include <algorithm>
#include <boost/multi_array.hpp>
#include <boost/shared_array.hpp>
#include <cstddef>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <iomanip>

#include "HungarianMethod.hpp"

namespace rfs{

  /**
   * \class OSPA
   * Optimal Sub-pattern Assignment metric
   * \f[ e_{\textrm{OSPA}} = \Bigg(\frac{1}{\left|\mathcal{M}_g\right|} \Bigg[ \min_{j \in \left\{ 1 \ldots \left|\mathcal{M}_g\right|\right\}} \sum_{i=1}^{\left|\mathcal{M}_k\right|} \min \left(d_{i,j},c\right)^p + c^p \Big| \left|\mathcal{M}_g\right| - \left|\mathcal{M}_k\right| \Big| \Bigg] \Bigg)^{\frac{1}{p}} \f]
   * \brief A class for computing the optimal sub-pattern assignment metric given 2 sets
   * \tparam SetElementType The class representing an element of the sets to be compared. The operator - must be defined in this class to return the distance between elements
   * \author Keith Leung
   */
  template<class SetElementType> 
  class OSPA{

  public:

    /** Cost matrix */
    typedef boost::multi_array<double, 2> CostMat;

    /** Cost matrix index */
    typedef CostMat::index Idx;

    /** Solution array */
    typedef boost::shared_array<int> SolnArr;

    /** Constructor 
     * \param set1 a std::vector containing the elements in set 1
     * \param set2 a std::vector containing the elements in set 2
     * \param cutoff value c
     * \param order value p
     */
    OSPA(std::vector<SetElementType> &set1,
	 std::vector<SetElementType> &set2,
	 double cutoff,
	 double order);

    /** Destructor */
    ~OSPA();

    /** Calculate the OSPA error
     * \param[out] e_dist distance error component (before raising to power p_ and averaging over n_) 
     * \param[out] e_card cardinality error component (before raising to power p_ and averaging over n_)
     * \param[in] report True to generate report to stdout
     * \return OSPA error
     */
    double calcError(double *e_dist = NULL, double *e_card = NULL, bool report = false);

    /** Report the optimal assignment using std::cout */
    void reportSoln();

    /** Get all the optimal assignments */
    SolnArr getOptAssignment();

    /** Get the cost matrix */
    CostMat getCostMat();

    /** Get an optimal assignment 
     *  \param[in] set1Idx index for set 1
     *  \param[out] cost The cost associated with the assignment. 
     *  \return the corresponding assignment from set 2, or -1 if not assigned or if set1Idx is invalid .
     */
    int getOptAssignment(unsigned int set1Idx, double *cost);

  protected:
    
    double c_; /**< cutoff threshold */
    double p_; /**< order of metric */
    CostMat C_; /**< Cost matrix */
    unsigned int n1_; /**< size of set 1 */
    unsigned int n2_; /**< size of set 2 */
    unsigned int n_; /**< size of cost matrix */
    SolnArr soln_; /**< Solution array */
    double cost_; /**< error */

  };

  // Implementation

  template<class SetElementType>
  OSPA<SetElementType>::OSPA(std::vector<SetElementType> &set1,
			     std::vector<SetElementType> &set2,
			     double cutoff,
			     double order) : c_(cutoff),
					     p_(order)
							      
  {

    // Allocate memory for cost matrix
    n1_ = set1.size();
    n2_ = set2.size();
    n_ = std::max(n1_, n2_);

    // Construct cost matrix
    C_.resize(boost::extents[n_][n_]); 
    for( Idx i = 0; i < n1_; i++ ){
      for( Idx j = 0; j < n2_; j++){ 
	C_[i][j] = fabs(set1[i] - set2[j]); 
	if(C_[i][j] > c_)
	  C_[i][j] = c_;
      }
    }
    if( n1_ < n_ ){ // n1 < n2 
      for( Idx i = n1_; i < n_; i++){
	for( Idx j = 0; j < n_; j++){
	  C_[i][j] = c_;
	}
      }
    }else if( n2_ < n_ ){ // n2 < n1
      for( Idx i = 0; i < n_; i++){
	for( Idx j = n2_; j < n_; j++){
	  C_[i][j] = c_;
	}
      }
    }

    /*for( Idx i = 0; i < n_; i++ ){
      for( Idx j = 0; j < n_; j++){ 
	std::cout << std::setw(10) << C_[i][j];
      }
      std::cout << std::endl;
      }*/

    // Use the Hungarian method to find optimal assignment
    rfs::HungarianMethod hm;
    soln_.reset( new int[n_] );
    double cost; // This cost is not the ospa cost because order p is not accounted for
    hm.run<CostMat>(C_, n_, soln_.get(), &cost, false);

  }

  template<class SetElementType>
  OSPA<SetElementType>::~OSPA(){}

  template<class SetElementType>
  double OSPA<SetElementType>::calcError(double *e_dist, double *e_card, bool report){
 
    cost_ = 0;
    if(e_dist != NULL)
      *e_dist = 0;
    if(e_card != NULL)
      *e_card = 0;
  
    for(Idx i = 0; i < n_; i++){
      unsigned int j = soln_[i];
      if( C_[i][j] == c_ && e_card != NULL){
	*e_card += C_[i][j];
      }else if(e_dist != NULL){
	*e_dist += C_[i][j];
      }
      cost_ += pow(C_[i][j], p_);
    }
    cost_ = pow(cost_ / n_, 1.0 / p_);

    return cost_;
  }
  
  template<class SetElementType>
  void OSPA<SetElementType>::reportSoln(){

    printf("Assignment Results:\n\n");

    for(Idx i = 0; i < n_; i++){
      Idx j = soln_[i];
      if(i >= n1_){
	printf("n/a -- %03ld [%lf]\n", j, C_[i][j]);
      }else if(j >= n2_){
	printf("%03ld -- n/a [%lf]\n", i, C_[i][j]);
      }else{
	printf("%03ld -- %03ld [%lf]\n", i, j, C_[i][j]);
      }
    }
    calcError();
    printf("\nError: %f\n", cost_);
    
  }

  template<class SetElementType>
  typename OSPA<SetElementType>::SolnArr OSPA<SetElementType>::getOptAssignment(){

    return soln_; 

  }

  template<class SetElementType>
  int OSPA<SetElementType>::getOptAssignment(unsigned int set1Idx, double *cost){

    if(set1Idx >= n1_)
      return -1;
    Idx j = soln_[set1Idx];
    if(j >= n2_){
      return -1;
    }else{
      *cost = C_[set1Idx][j];
      return j;
    }
    
  }
  
  template<class SetElementType>
  typename OSPA<SetElementType>::CostMat OSPA<SetElementType>::getCostMat(){
    return C_;
  }

} // namespace rfs

#endif
