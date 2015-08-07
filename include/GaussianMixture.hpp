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

#ifndef GAUSSIAN_MIXTURE_HPP
#define GAUSSIAN_MIXTURE_HPP

#include <algorithm>
#include <iostream>
#include "Landmark.hpp"
#include <vector>

namespace rfs{

/** 
 * \class GaussianMixture
 * \brief A class for mixture of Gaussians
 *
 * This class represents a mixture of weighted Gaussians
 * It is designed to work with the RBPHDFilter
 *
 * \author Keith Leung
 */

template< class Landmark >
class GaussianMixture
{
public:

  typedef Landmark TLandmark;
  typedef Landmark* pLandmark;

  /** \brief A data structure representing a weighted Gaussian distribution in GaussianMixture */ 
  struct Gaussian{
    pLandmark landmark; /**< pointer to a Landmark, which holds the mean and covariance */
    double weight; /**< weight of Gaussian */
    double weight_prev; /**< previous weight of Gaussian (used in computations by the RBPHDFilter) */
  };
  
  /** Default constructor */
  GaussianMixture();

  /** Copy constructor */
  GaussianMixture(const GaussianMixture& other);
  
  /** Destructor */
  ~GaussianMixture();

  /** 
   * Copy data from this GaussianMixture to another GaussianMixture.
   * Memory is allocated for creating copies of Gaussians for the other GaussianMixture.
   * \param[in,out] other the other Gaussian mixture to which data is copied to
   */
  void copyTo( GaussianMixture *other);

  /** 
   * Add a Gaussian to this GaussianMixture
   * \param[in] p pointer to the Landmark which holds the mean and covariance
   * \param[in] w weight of the new Gaussian
   * \param[in] allocateMem if false, the function assumes memory for Landmark has already been allocated
   * and will not go out of scope or get deleted, other than by the current instantiation of GaussianMixture.
   * If true, memory is allocated for a new Landmark and data from p is copied to it.
   * \return number of Gaussians in the mixture
   */
  unsigned int addGaussian( pLandmark p, double w = 1, bool allocateMem = false);

  /** 
   * Remove a Gaussian from the mixture 
   * \param[in] idx index number of the Gaussian to remove
   * \return number of Gaussians in the mixture
   */
  unsigned int removeGaussian( unsigned int idx );

  /**
   * Get the number of Gaussians in the mixture
   * \return count 
   */
  unsigned int getGaussianCount();

  /** 
   * Get the sum of all Gaussian weights 
   * \return sum
   */
  double getWeightSum();

  /**
   * Set the weight of a Gaussian indicated by the given index
   * \param[in] idx index
   * \param[in] w weight 
   */ 
  void setWeight( unsigned int idx, double w );

  /**
   * Get the weight of a Gaussian indicated by the given index
   * \param[in] idx index
   * \return weight 
   */ 
  double getWeight( unsigned int idx );

  /**
   * Get the parameters of a Gaussian (stored as a Landmark)
   * \param[in] idx index
   * \return pointer to a Landmark, or NULL if the Gaussian does not exist.
   */
  pLandmark getGaussian( unsigned int idx );
  
  /**
   * Get the parameters of a Gaussian (stored as a Landmark)
   * \param[in] idx index
   * \param[out] p pointer to a Landmark that holds the Gaussian parameters
   */
  void getGaussian( unsigned int idx, pLandmark &p);

  /**
   * Get the state and weight of a Gaussian
   * \param[in] idx index
   * \param[out] p pointer to a Landmark which holds the mean and covariance
   * \param[out] w the weight 
   */
  void getGaussian( unsigned int idx, pLandmark &p, double &w);

  /**
   * Get the parameters and weight of a Gaussian
   * \param[in] idx index
   * \param[out] p pointer to a Landmark which holds the mean and covariance
   * \param[out] w overwritten by the weight 
   * \param[out] w_prev overwritten by the previous weight (used in parts of the RBPHDFilter)
   */
  void getGaussian( unsigned int idx, pLandmark &p, double &w, double &w_prev);

  /**
   * Update an Gaussian in the mixture. 
   * \param[in] idx index of the Gaussian to update
   * \param[in] lm Landmark object with the updated mean and covariance
   * \param[in] w weight of the updated Gaussian. No change are made to the existing weight if this is negative.
   * \return true if idx is valid and Gaussian is updated
   */
  bool updateGaussian( unsigned int idx, Landmark &lm, double w = -1);

  /**
   * Check all Gaussians in the mixture and merge those that are within a certain Mahalanobis distance of each other
   * \param[in] t distance threshold (the default value squares to 0.1)
   * \param[in] f_inflation merged Gaussian covariance inflation factor (default value causes no inflation)
   * \return number of merging operations
   */
  unsigned int merge(const double t = 0.31622776601, const double f_inflation = 1.0);

  /**
   * Prune the Gaussian mixture to remove Gaussians with weights that are
   * less than a given threshold. 
   * \note Gaussians may have difference indices after using this function due to sorting performed on gList_
   * \param[in] t weight threshold, below which Gaussians are removed.
   * \return number of Gaussians removed
   */
  unsigned int prune( const double t );

  /**					       
   * Sort the Gaussian mixture container gList_ from highest to lowest Gaussian weight.
   */
  void sortByWeight();

protected:

  int n_; /**< number of Gaussians in this GaussianMixture */
  std::vector<Gaussian> gList_; /**< container for Gaussians */ 
  bool isSorted_; /**< a flag to prevent unecessary sorting */
  
  /**
   * Comparison function used by sortByWeight() for sorting the Gaussian container gList_.
   */
  static bool weightCompare(Gaussian a, Gaussian b);

  /** 
   * Add a Gaussian to this GaussianMixture, by overwrite an existing spot in the Gaussian container.
   * \param[in] idx element in the container gList_ to overwrite.
   * \param[in] p pointer to the Landmark to add, which holds the mean and covariance.
   * \param[in] w weight of the new Gaussian.
   * \param[in] allocateMem if false, this assumes memory for Landmark has already been allocated
   * and will not go out of scope or get deleted, other than by the current instantiation of GaussianMixture.
   * If true, memory is allocated for a new Landmark and data from p is copied to it.
   * \return number of Gaussians in the mixture
   */  
  unsigned int addGaussian( unsigned int idx, pLandmark p, double w = 1, bool allocateMem = false);

  /**
   * Merge two Guassians if the second is within a Mahalanobis distance of the first. 
   * If merging occurs, the resulting Gaussian will overwrite the first Gaussian,
   * and second one will be removed from the Gaussian mixture.
   * \param[in] idx1 index of the first Gaussian
   * \param[in] idx2 index of the second Gaussian
   * \param[in] t distance threshold (the default value squares to 0.1)
   * \param[in] f_inflation merged Gaussian covariance inflation factor
   * \return true if merging is successful
   */
  bool merge(unsigned int idx1, unsigned int idx2, 
	     const double t = 0.3162277660, const double f_inflation = 1.0);

};

////////// Implementation //////////

template< class Landmark >
GaussianMixture<Landmark>::GaussianMixture(){
  n_ = 0;
  isSorted_ = false;
}

template< class Landmark >
GaussianMixture<Landmark>::GaussianMixture(const GaussianMixture& other):
n_(other.n_), isSorted_(other.isSorted_)
{
  gList_ = other.gList_;
  for(int i = 0; i < gList_.size(); i++){
    Landmark* lmkCopy = new Landmark;
    *lmkCopy = *(other.gList_[i].landmark);
    gList_[i].landmark = lmkCopy; 
  }
}

template< class Landmark >
GaussianMixture<Landmark>::~GaussianMixture(){

  for( int i = 0; i < gList_.size(); i++){
    removeGaussian(i);
  }

}

template< class Landmark >
void GaussianMixture<Landmark>::copyTo( GaussianMixture *other){
  other->n_ = n_;
  other->gList_ = gList_;
  for(int i = 0; i < gList_.size(); i++){
    Landmark* lmkCopy = new Landmark;
    *lmkCopy = *(gList_[i].landmark);
    other->gList_[i].landmark = lmkCopy; 
  }
  other->isSorted_ = isSorted_;
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::addGaussian( pLandmark p, double w, 
						     bool allocateMem){
  isSorted_ = false;

  Landmark* pNew = NULL;
  if(allocateMem){
    pNew = new Landmark;
    *pNew = *p;
  }else{
    pNew = p;
  }

  Gaussian g = {pNew, w, 0};
  gList_.push_back(g);
  n_++;
  return n_;
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::addGaussian( unsigned int idx, pLandmark p, double w, 
						     bool allocateMem){
  isSorted_ = false;

  Landmark* pNew = NULL;
  if(allocateMem){
    pNew = new Landmark;
    *pNew = *p;
  }else{
    pNew = p;
  }

  if (gList_[idx].landmark != NULL){
    removeGaussian( idx );
  }

  gList_[idx].landmark = pNew;
  gList_[idx].weight = w;
  gList_[idx].weight_prev = 0;
  n_++;
  return n_;
}

template< class Landmark > 
unsigned int GaussianMixture<Landmark>::removeGaussian( unsigned int idx ){

  isSorted_ = false;
  if(gList_[idx].landmark != NULL){
    delete gList_[idx].landmark;
    gList_[idx].landmark = NULL;
    gList_[idx].weight = 0;
    gList_[idx].weight_prev = 0;
    n_--;
  }
  return n_;
}


template< class Landmark >
unsigned int GaussianMixture<Landmark>::getGaussianCount(){
  return n_;
}

template< class Landmark >
double GaussianMixture<Landmark>::getWeightSum(){
  double sum = 0;
  for(int i = 0; i < n_; i++){
    sum += gList_[i].weight;
  }
  return sum;
}

template< class Landmark >
void GaussianMixture<Landmark>::setWeight( unsigned int idx, double w ){
  isSorted_ = false;
  gList_[idx].weight_prev = gList_[idx].weight;
  gList_[idx].weight = w;
}

template< class Landmark >
double GaussianMixture<Landmark>::getWeight( unsigned int idx){
  return gList_[idx].weight;
}

template< class Landmark >
typename GaussianMixture<Landmark>::pLandmark GaussianMixture<Landmark>::getGaussian( unsigned int idx ){
  return (gList_[idx].landmark);
}

template< class Landmark >
void GaussianMixture<Landmark>::getGaussian( unsigned int idx, pLandmark &p){
  p = gList_[idx].landmark; 
}

template< class Landmark >
void GaussianMixture<Landmark>::getGaussian( unsigned int idx, pLandmark &p, double &w){
    p = gList_[idx].landmark;
    w = gList_[idx].weight;
}

template< class Landmark >
void GaussianMixture<Landmark>::getGaussian( unsigned int idx, pLandmark &p, double &w, double &w_prev){
    p = gList_[idx].landmark;
    w = gList_[idx].weight;
    w_prev = gList_[idx].weight_prev;
}


template< class Landmark >
bool GaussianMixture<Landmark>::updateGaussian( unsigned int idx, Landmark &lm, double w){
  
  isSorted_ = false;

  if( idx > gList_.size() || gList_[idx].landmark == NULL) 
    return false;

  if ( w < 0 )
    w = getWeight( idx );

  *(gList_[idx].landmark) = lm;
  gList_[idx].weight_prev = gList_[idx].weight;
  gList_[idx].weight = w;

  return true;

  }

template< class Landmark >
unsigned int GaussianMixture<Landmark>::merge(const double t, 
					      const double f_inflation){
  isSorted_ = false;
  unsigned int nMerged = 0;
  int nGaussians = gList_.size();

  for(unsigned int i = 0; i < nGaussians; i++){
   
    if(gList_[i].landmark == NULL)
      continue;

    for(unsigned int j = i+1; j < nGaussians; j++){

      if( merge(i, j, t, f_inflation) ){
	nMerged++;
      } 
	
    }
    
  }
  return nMerged;
}


template< class Landmark >
bool GaussianMixture<Landmark>::merge(unsigned int idx1, unsigned int idx2,
				      const double t, const double f_inflation){
  isSorted_ = false;

  if (gList_[idx1].landmark == NULL ||
      gList_[idx2].landmark == NULL ){
    return false;
  }

  double w_m, w_1, w_2; 
  double d_mahalanobis_1, d_mahalanobis_2; 

  w_1 = gList_[idx1].weight;
  w_2 = gList_[idx2].weight;
  
  double t2 = t * t;
  d_mahalanobis_1 = gList_[idx1].landmark->mahalanobisDist2( *(gList_[idx2].landmark) );
  if( d_mahalanobis_1 > t2 ){
    d_mahalanobis_2 = gList_[idx2].landmark->mahalanobisDist2( *(gList_[idx1].landmark) );
    if( d_mahalanobis_2 > t2 ){
      return false;
    }
  }

  w_m = w_1 + w_2;
  
  if( w_m == 0 )
    return false;

  typename Landmark::Vec x_1, x_2, x_m, d_1, d_2;
  typename Landmark::Mat S_1, S_2, S_m;

  gList_[idx1].landmark->get(x_1, S_1);
  gList_[idx2].landmark->get(x_2, S_2);

  x_m = (x_1 * w_1 + x_2 * w_2) / w_m;

  d_1 = x_m - x_1;
  d_2 = x_m - x_2;
  S_m = ( w_1 * ( S_1 + f_inflation * d_1 * d_1.transpose() ) +
          w_2 * ( S_2 + f_inflation * d_2 * d_2.transpose() ) ) / w_m;


  //std::cout << "\n" << x_1 << "\n" << S_1 << "\n\n";
  //std::cout << "\n" << x_2 << "\n" << S_2 << "\n\n";
  //std::cout << "\n" << x_m << "\n" << S_m << "\n\n";

  gList_[idx1].landmark->set(x_m, S_m);
  gList_[idx1].weight = w_m;
  gList_[idx1].weight_prev = 0;

  removeGaussian( idx2 ); // this also takes care of updating Gaussian count

  return true;
  
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::prune( const double t ){

  isSorted_ = false;

  unsigned int nPruned = 0;
  if( gList_.size() < 1 ){
    return nPruned;
  }

  sortByWeight(); // Sort from greatest to smallest weight
  
  // Binary search for the Gaussian with weight closest to and greater than t
  unsigned int min_idx = 0;
  unsigned int max_idx = gList_.size() - 1;
  unsigned int idx = (unsigned int)( (max_idx + min_idx) / 2 );
  unsigned int idx_old = idx + 1;
  double w = gList_[idx].weight;

  while(idx != idx_old){
    if( w <= t ){ 
      max_idx = idx;
    }else if ( w > t ){
      min_idx = idx;
    }
    idx_old = idx;
    idx = (unsigned int)( (max_idx + min_idx) / 2 );
    w = gList_[idx].weight;
  }
  while( w >= t ){
    idx++;
    if(idx >= gList_.size())
      break;
    w = gList_[idx].weight;
  }
  idx_old = idx; 
  while( idx < gList_.size() ){
    removeGaussian(idx); // this already takes care updating the Gaussian count
    idx++;
    nPruned++;
  }
  gList_.resize( gList_.size() - nPruned); 

  return nPruned;
}

template< class Landmark >
void GaussianMixture<Landmark>::sortByWeight(){
  if( !isSorted_ ){
    std::sort( gList_.begin(), gList_.end(), weightCompare );
    isSorted_ = true;
  }
}

template< class Landmark >
bool GaussianMixture<Landmark>::weightCompare(Gaussian a, Gaussian b){
  return a.weight > b.weight;
}

}

#endif
