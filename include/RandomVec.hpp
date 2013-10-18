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

// State class
// Keith Leung 2013

#ifndef STATE_HPP
#define STATE_HPP


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
//#include <Eigen/StdVector>
#include <iostream>
#include <stdio.h>

double const PI = acos(-1);

/**
 * \class RandomVec
 * Representation of a Gaussian random vector, with mean and covariance. A time variable is also included for time-stamping.
 * \brief An abstract base class for deriving pose and measurement classes
 * \tparam VecType An Eigen vector of dimension n
 * \tparam MatType An Eigen matrix or dimension n x n
 * \author Keith Leung
 */
template<class VecType, class MatType>
class RandomVec
{

public:

  typedef VecType Vec;
  typedef MatType Mat;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** Default constructor */
  RandomVec() : 
    isValid_Sx_L_(false), 
    isValid_Sx_inv_(false),
    isValid_Sx_det_(false),
    gen_(NULL)
  {

    if ( !dimCheck() ){
      exit(-1);
    }
    x_.setZero();
    Sx_.setZero();
    t_ = -1;

  }

  /** 
   * Constructor 
   * \param[in] x vector
   * \param[in] Sx covariance
   * \param[in] t time
   */
  RandomVec(VecType x, MatType Sx, double t = -1) : 
    isValid_Sx_L_(false), 
    isValid_Sx_inv_(false),
    isValid_Sx_det_(false),
    gen_(NULL)
  {
    if ( !dimCheck() ){
      exit(-1);
    }
    set(x);
    setCov(Sx);
    t_ = t;

  }

  /** 
   * Constructor 
   * \param[in] x vector
   * \param[in] t time
   */
  RandomVec(VecType x, double t = -1) : 
    isValid_Sx_L_(false), 
    isValid_Sx_inv_(false),
    isValid_Sx_det_(false),
    gen_(NULL)
  {
    if ( !dimCheck() ){
      exit(-1);
    }
    set(x);
    Sx_.setZero();
    t_ = t;

  }

  
  /**
   * Copy constructor
   * \param[in] other the RandomVec being copied from
   */
  RandomVec( const RandomVec& other ):

    x_( other.x_ ), 
    nDim_( other.nDim_ ), 
    Sx_( other.Sx_ ), 
    Sx_inv_( other.Sx_inv_ ),
    isValid_Sx_inv_( other.isValid_Sx_inv_),
    Sx_det_( other.Sx_det_ ),
    isValid_Sx_det_( other.isValid_Sx_det_),
    Sx_L_( other.Sx_L_ ),
    isValid_Sx_L_( other.isValid_Sx_L_), 
    t_(other.t_) 
  {
    gen_ = NULL;
  }

  /**
   * Assignment operator
   * \param[in] rhs the right-hand-side from which data is copied
   */
  RandomVec& operator=( const RandomVec& rhs ){
    
    x_ = rhs.x_;
    nDim_ = rhs.nDim_;
    Sx_ = rhs.Sx_;
    Sx_inv_ = rhs.Sx_inv_;
    isValid_Sx_inv_ = rhs.isValid_Sx_inv_;
    Sx_det_ = rhs.Sx_det_;
    isValid_Sx_det_ = rhs.isValid_Sx_det_;
    Sx_L_ = rhs.Sx_L_;
    isValid_Sx_L_ = rhs.isValid_Sx_L_;
    t_ = rhs.t_;

    gen_ = NULL;
  }

  /** Default destructor */
  ~RandomVec(){
    if( gen_ != NULL )
      delete gen_;
  };

  /** 
   * Set the vector
   * \param[in] x vector to be set
   */
  void set( VecType &x ){x_ = x;}

  /** 
   * Set the covariance for uncertainty
   * \param[in] Sx uncertainty to be set
   */
  void setCov( MatType &Sx){
    Sx_ = Sx;
    isValid_Sx_L_ = false; 
    isValid_Sx_inv_ = false;
    isValid_Sx_det_ =false;
  }

  /**
   * Set the time
   * \param[in] t time
   */
  void setTime( double t ){
    t_ = t;
  }

  /** 
   * Set the vector with a covariance matrix
   * \param[in] x vector to be set
   * \param[in] Sx covariance to be set
   */
  void set( VecType &x, MatType &Sx){
    set(x);
    setCov(Sx);
  }

  /** 
   * Set the vector with a time
   * \param[in] x vector to be set
   * \param[in] t time
   */
  void set( VecType &x, double t){
    set(x);
    t_ = t;
  }

  /** 
   * Set the vector with a covariance matrix, and time
   * \param[in] x vector to be set
   * \param[in] Sx covariance to be set
   * \param[in] t time
   */
  void set( VecType &x, MatType &Sx, double t){
    set(x);
    setCov(Sx);
    t_ = t;
  }


  /** 
   * Get the vector
   * \param[out] x vector
   */
  void get( VecType &x ){x = x_;}

  /** 
   * Get the covariance matrix
   * \param[out] Sx uncertainty 
   */
  void getCov( MatType &Sx){
    Sx = Sx_;
  }

  /**
   * Get the lower triangular part of the Cholesky decomposition 
   * on the covariance matrx Sx_ 
   * \param[out] Sx_Chol_L the lower triangular part of the Choloesky decomposition
   */
  void getCovCholeskyDecompLower( MatType &Sx_Chol_L){
    if(!isValid_Sx_L_){
      Eigen::LLT<MatType> cholesky( Sx_ );
      Sx_L_ = cholesky.matrixL();
      isValid_Sx_L_ = true;
    }
    Sx_Chol_L = Sx_L_;
  }

  /** 
   * Get the invserve covariance matrix 
   * \param[out] Sx_inv inverse covariance
   */ 
  void getCovInv( MatType &Sx_inv){
    if(!isValid_Sx_inv_){
      Sx_inv_ = Sx_.inverse(); 
      isValid_Sx_inv_ = true;
    }
    Sx_inv = Sx_inv_;
  }

  /** 
   * Get the determinant of the covariance
   * \return determinant
   */
  double getCovDet(){
    if(!isValid_Sx_det_){
      Sx_det_ = Sx_.determinant();
      isValid_Sx_det_ = true;
    }
    return Sx_det_;
  }

  /** 
   * Get the vector and covariance matrix
   * \param[out] x vector
   * \param[out] Sx uncertainty
   */
  void get( VecType &x, MatType &Sx){
    get(x);
    getCov(Sx);
  }

  /** 
   * Get the vector and time
   * \param[out] x vector
   * \param[out] t time
   */
  void get( VecType &x, double &t){
    get(x);
    t = t_;
  }

  /** 
   * Get the vector, covariance matrix, and time
   * \param[out] x vector
   * \param[out] Sx uncertainty
   * \param[out] t time
   */
  void get( VecType &x, MatType &Sx, double &t){
    get(x);
    getCov(Sx);
    t = t_;
  }

  /** 
   * Get an element of the vector
   * \param[in] n element index
   * \return element n
   */
  double get( int n ){ return x_(n);}

  /**
   * Get the time
   * \return time
   */
  double getTime(){
    return t_;
  }

  /** 
   * Get the dimension
   * \return dimension
   */ 
  unsigned int getNDim(){ return nDim_; }

  /**
   * Calculate the squared Mahalanobis distance to another random vector of the same type
   * \param[in] to the RandomVec containing the vector we are measuring to 
   */
  double mahalanobisDist2( RandomVec<VecType, MatType> &to ){
    if(!isValid_Sx_inv_){
      Sx_inv_ = Sx_.inverse(); 
      isValid_Sx_inv_ = true;
    }
    e_ = to.x_ - x_;
    return (e_.transpose() * Sx_inv_ * e_);
  }

  /**
   * Calculate the squared Mahalanobis distance to another random vector of the same type
   * \param[in] to_x the vector we are measuring to 
   */
  double mahalanobisDist2( typename RandomVec<VecType, MatType>::Vec &to_x ){
    if(!isValid_Sx_inv_){
      Sx_inv_ = Sx_.inverse();
      isValid_Sx_inv_ = true;
    }
    e_ = to_x - x_;
    return (e_.transpose() * Sx_inv_ * e_);
  }

  /**
   * Calculate the Gaussian likelihood of a given evaluation point
   * \param[in] x_eval the evaluation point
   * \param[out] if not NULL, the pointed to variable will be overwritten by the 
   * squared mahalanobis distance used to calculate the likelihood
   */ 
  double evalGaussianLikelihood( RandomVec<VecType, MatType> &x_eval,
				 double* mDist2 = NULL){
    if(!isValid_Sx_det_){
      Sx_det_ = Sx_.determinant();
      isValid_Sx_det_ = true;
    }
    double md2 = mahalanobisDist2( x_eval );
    double l = ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim_ ) * Sx_det_ ) );
    if( l != l) //If md2 is very large, l will become NAN;
      l = 0;
    if(mDist2 != NULL)
      *mDist2 = md2;
    return l;
  }

  /**
   * Calculate likelihood
   * \param[in] x_eval the evaluation point
   * \param[out] if not NULL, the pointed to variable will be overwritten by the 
   * squared mahalanobis distance used to calculate the likelihood
   */ 
  double evalGaussianLikelihood( typename RandomVec<VecType, MatType>::Vec &x_eval,
				 double* mDist2 = NULL){
    if(!isValid_Sx_det_){
      Sx_det_ = Sx_.determinant();
      isValid_Sx_det_ = true;
    }
    double md2 = mahalanobisDist2( x_eval );
    double l = ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim_ ) * Sx_det_ ) );
    if( l != l) //If md2 is very large, l will become NAN;
      l = 0;
    if(mDist2 != NULL)
      *mDist2 = md2;
    return l;
  }

  /** 
   * Sample this random vector
   * \param[out] s_sample The sampled vector, with time and covariance copied from this random vector
   */
  void sample( RandomVec<VecType, MatType> &s_sample ){
    
    VecType x_sample, indep_noise;

    if(!isValid_Sx_L_){
      Eigen::LLT<MatType> cholesky( Sx_ );
      Sx_L_ = cholesky.matrixL();
      isValid_Sx_L_ = true;
    }

    if(gen_ == NULL){
      gen_ = new boost::variate_generator< boost::mt19937, 
					   boost::normal_distribution<double> >
	(boost::mt19937(rand()), boost::normal_distribution<double>());
    }
    
    int n = Sx_L_.cols();
    for(int i = 0; i < n; i++){
      indep_noise(i) = (*gen_)();
    }
    x_sample = x_ + Sx_L_ * indep_noise;
    s_sample.set( x_sample, Sx_, t_ );

  }


  /** 
   * Sample this random vector and write the result over the mean x_;
   */
  void sample(){
    
    VecType x_sample, indep_noise;

    if(!isValid_Sx_L_){
      Eigen::LLT<MatType> cholesky( Sx_ );
      Sx_L_ = cholesky.matrixL();
      isValid_Sx_L_ = true;
    }

    if(gen_ == NULL){
      gen_ = new boost::variate_generator< boost::mt19937, 
					   boost::normal_distribution<double> >
	(boost::mt19937(rand()), boost::normal_distribution<double>());
    }
    
    int n = Sx_L_.cols();
    for(int i = 0; i < n; i++){
      indep_noise(i) = (*gen_)();
    }
    x_ += Sx_L_ * indep_noise;

  }

private:

  VecType x_; /**< State */
  unsigned int nDim_; /**< Number of dimensions */
  MatType Sx_; /**< Covariance */
  MatType Sx_inv_; /**< Inverse covariance */
  bool isValid_Sx_inv_; /**< Inverse covariance is up to date */
  double Sx_det_; /**< Determinant of Sx_ */
  bool isValid_Sx_det_; /**< Determinant of Sx_ is up to date */
  MatType Sx_L_; /**< Lower triangular part of Cholesky decomposition on Sx_ */
  bool isValid_Sx_L_; /**< Lower triangular part of Cholesky decomposition on Sx_ is up to date */
  double t_; /**< time */

  VecType e_; /**< temporary */

  boost::variate_generator< boost::mt19937, 
			    boost::normal_distribution<double> >* gen_;/**< normal distribution random number generator */ 

  /** Dimensionality check during initialization */
  bool dimCheck(){

    if( Sx_.rows() != Sx_.cols() ){
      std::cerr << "Error: MatType must be a square matrix \n";
      return false;
    }
    if( Sx_.rows() != x_.size() ){
      std::cerr << "Error: VecType and MatType dimension mismatch \n";
      return false;
    }
    nDim_ = x_.size();
    return true;
  }

};

#endif
