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

#include <assert.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/LU>
//#include <Eigen/StdVector>
#include <iostream>
#include <stdio.h>
#include "TimeStamp.hpp"


namespace rfs
{


  double const PI = acos(-1);

  /**
   * \class RandomVec
   * Representation of a Gaussian random vector, with mean and covariance. A time variable is also included for time-stamping.
   * \brief An abstract base class for deriving pose and measurement classes
   * \tparam nDim Dimension of the vector
   * \author Keith Leung
   */
  template<unsigned int nDim = 1>
  class RandomVec
  {

  public:

    typedef ::Eigen::Matrix<double, nDim, 1> Vec;
    typedef ::Eigen::Matrix<double, nDim, nDim> Mat;
    typedef Mat Cov;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    /** \brief Default constructor */
    RandomVec() : 
      isValid_Sx_L_(false), 
      isValid_Sx_inv_(false),
      isValid_Sx_det_(false)
    {
      dimCheck();

      x_.setZero();
      Sx_.setZero();

    }

    /** 
     * \brief Constructor 
     * \param[in] x vector
     * \param[in] Sx covariance
     * \param[in] t time
     */
    RandomVec(const Vec &x, const Mat &Sx, const TimeStamp &t = TimeStamp() ) :
      isValid_Sx_L_(false), 
      isValid_Sx_inv_(false),
      isValid_Sx_det_(false)
    {
      dimCheck();
      set(x);
      setCov(Sx);
      t_ = t;

    }

    /** 
     * \brief Constructor 
     * \param[in] x vector
     * \param[in] SxVec Array with diagonal entries of covariance
     * \param[in] t time
     */
    RandomVec(const Vec &x, double const * const &SxVec, const TimeStamp &t = TimeStamp() ) :
      isValid_Sx_L_(false), 
      isValid_Sx_inv_(false),
      isValid_Sx_det_(false)
    {
      dimCheck();
      set(x);
      Vec SxVecTmp;
      for(int i = 0; i < nDim; i++){
	SxVecTmp(i) = SxVec[i]; 
      }
      setCov(SxVecTmp.asDiagonal());
      t_ = t;

    }

    /** 
     * Constructor 
     * \param[in] x vector
     * \param[in] t time
     */
    RandomVec(const Vec &x, const TimeStamp t = TimeStamp()) :
      isValid_Sx_L_(false), 
      isValid_Sx_inv_(false),
      isValid_Sx_det_(false)
    {
      dimCheck();
      set(x);
      Sx_.setZero();
      t_ = t;

    }

    /** 
     * Constructor 
     * \param[in] t time
     */
    RandomVec(const TimeStamp &t) :
      isValid_Sx_L_(false), 
      isValid_Sx_inv_(false),
      isValid_Sx_det_(false)
    {
      dimCheck();
      x_.setZero();
      Sx_.setZero();
      t_ = t;

    }
  
    /**
     * Copy constructor
     * \param[in] other the RandomVec being copied from
     */
    RandomVec( const RandomVec& other ):
    x_( other.x_ ), 
    Sx_( other.Sx_ ), 
    Sx_inv_( other.Sx_inv_ ),
    isValid_Sx_inv_( other.isValid_Sx_inv_),
    Sx_det_( other.Sx_det_ ),
    isValid_Sx_det_( other.isValid_Sx_det_),
    Sx_L_( other.Sx_L_ ),
    isValid_Sx_L_( other.isValid_Sx_L_), 
    t_(other.t_) 
    {}

    /**
     * Assignment operator
     * \param[in] rhs the right-hand-side from which data is copied
     */
    RandomVec& operator=( const RandomVec& rhs ){
    
      x_ = rhs.x_;
      Sx_ = rhs.Sx_;
      Sx_inv_ = rhs.Sx_inv_;
      isValid_Sx_inv_ = rhs.isValid_Sx_inv_;
      Sx_det_ = rhs.Sx_det_;
      isValid_Sx_det_ = rhs.isValid_Sx_det_;
      Sx_L_ = rhs.Sx_L_;
      isValid_Sx_L_ = rhs.isValid_Sx_L_;
      t_ = rhs.t_;

      return *this;
    }

    /** Default destructor */
    ~RandomVec(){};

    /** 
     * [] Operator for looking up the value of an element of x			       
     * \return the value of element n of vector x
     */
    double& operator[] (const int n){ 
      assert(n >= 0 && n < nDim);
      return x_(n);
    }

    /** 
     * Set the vector
     * \param[in] x vector to be set
     */
    void set( const Vec &x ){x_ = x;}

    /** 
     * Set the covariance for uncertainty
     * \param[in] Sx uncertainty to be set
     */
    void setCov( const Mat &Sx){
      Sx_ = Sx;
      isValid_Sx_L_ = false; 
      isValid_Sx_inv_ = false;
      isValid_Sx_det_ =false;
    }

    /**
     * Set the time
     * \param[in] t time
     */
    void setTime( const TimeStamp &t ){
      t_ = t;
    }

    /** 
     * Set the vector with a covariance matrix
     * \param[in] x vector to be set
     * \param[in] Sx covariance to be set
     */
    void set( const Vec &x, const Mat &Sx){
      set(x);
      setCov(Sx);
    }

    /** 
     * Set the vector with a time
     * \param[in] x vector to be set
     * \param[in] t time
     */
    void set( const Vec &x, const TimeStamp &t){
      set(x);
      t_ = t;
    }

    /** 
     * Set the vector with a covariance matrix, and time
     * \param[in] x vector to be set
     * \param[in] Sx covariance to be set
     * \param[in] t time
     */
    void set( const Vec &x, const Mat &Sx, const TimeStamp &t){
      set(x);
      setCov(Sx);
      t_ = t;
    }

    /**
     * Get the vector
     * \return x vector
     */
    Vec get() const { return x_;}

    /** 
     * Get the vector
     * \param[out] x vector
     */
    void get( Vec &x ) const {x = x_;}

    /**
     * Get the covariance matrix
     * \return Sx covariance representing the uncertainty
     */
    Mat getCov() const { return Sx_; }

    /** 
     * Get the covariance matrix
     * \param[out] Sx uncertainty representing the uncertainty 
     */
    void getCov( Mat &Sx) const {
      Sx = Sx_;
    }

    /**
     * Get the lower triangular part of the Cholesky decomposition 
     * on the covariance matrx Sx_ 
     * \param[out] Sx_Chol_L the lower triangular part of the Choloesky decomposition
     */
    void getCovCholeskyDecompLower( Mat &Sx_Chol_L){
      if(!isValid_Sx_L_){
	::Eigen::LLT<Mat> cholesky( Sx_ );
	Sx_L_ = cholesky.matrixL();
	isValid_Sx_L_ = true;
      }
      Sx_Chol_L = Sx_L_;
    }

    /** 
     * Get the invserve covariance matrix 
     * \param[out] Sx_inv inverse covariance
     */ 
    void getCovInv( Mat &Sx_inv){
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
    void get( Vec &x, Mat &Sx) const {
      get(x);
      getCov(Sx);
    }

    /** 
     * Get the vector and time
     * \param[out] x vector
     * \param[out] t time
     */
    void get( Vec &x, TimeStamp &t) const {
      get(x);
      t = t_;
    }

    /** 
     * Get the vector, covariance matrix, and time
     * \param[out] x vector
     * \param[out] Sx uncertainty
     * \param[out] t time
     */
    void get( Vec &x, Mat &Sx, TimeStamp &t) const {
      get(x);
      getCov(Sx);
      t = t_;
    }

    /** 
     * Get an element of the vector
     * \param[in] n element index
     * \return element n
     */
    double get( const int n ) const { return x_(n);}

    /**
     * Get the time
     * \return time
     */
    TimeStamp getTime() const {
      return t_;
    }

    /** 
     * Get the dimension
     * \return dimension
     */ 
    unsigned int getNDim() const { return nDim; }

    /**
     * Calculate the squared Mahalanobis distance to another random vector of the same type
     * \param[in] to the RandomVec containing the vector we are measuring to 
     */
    double mahalanobisDist2(const RandomVec<nDim> &to ){
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
    double mahalanobisDist2(const typename RandomVec<nDim>::Vec &to_x ){
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
     * \param[out] mDist2 if not NULL, the pointed to variable will be overwritten by the 
     * squared mahalanobis distance used to calculate the likelihood
     */ 
    double evalGaussianLikelihood(const RandomVec<nDim> &x_eval,
				  double* mDist2 = NULL){
      if(!isValid_Sx_det_){
	Sx_det_ = Sx_.determinant();
	gaussian_pdf_factor_ = sqrt( pow( 2*PI, nDim ) * Sx_det_ );
	isValid_Sx_det_ = true;
      }
      double md2 = mahalanobisDist2( x_eval );
      double l = ( exp(-0.5 * md2 ) / gaussian_pdf_factor_ );
      if( l != l) //If md2 is very large, l will become NAN;
	l = 0;
      if(mDist2 != NULL)
	*mDist2 = md2;
      return l;
    }

    /**
     * Calculate likelihood
     * \param[in] x_eval the evaluation point
     * \param[out] mDist2 if not NULL, the pointed to variable will be overwritten by the 
     * squared mahalanobis distance used to calculate the likelihood
     */ 
    double evalGaussianLikelihood(const typename RandomVec<nDim>::Vec &x_eval,
				  double* mDist2 = NULL){
      if(!isValid_Sx_det_){
	Sx_det_ = Sx_.determinant();
	gaussian_pdf_factor_ = sqrt( pow( 2*PI, nDim ) * Sx_det_ );
	isValid_Sx_det_ = true;
      }
      double md2 = mahalanobisDist2( x_eval );
      double l = ( exp(-0.5 * md2 ) / gaussian_pdf_factor_ );
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
    void sample( RandomVec<nDim> &s_sample ){
    
      Vec x_sample, indep_noise;

      if(!isValid_Sx_L_){
	::Eigen::LLT<Mat> cholesky( Sx_ );
	Sx_L_ = cholesky.matrixL();
	isValid_Sx_L_ = true;
      }

      int n = Sx_L_.cols();
      for(int i = 0; i < n; i++){
	indep_noise(i) = genGaussian_();
      }
      x_sample = x_ + Sx_L_ * indep_noise;
      s_sample.set( x_sample, Sx_, t_ );

    }


    /** 
     * Sample this random vector and write the result over the mean x_;
     */
    void sample(){
    
      Vec x_sample, indep_noise;

      if(!isValid_Sx_L_){
	::Eigen::LLT<Mat> cholesky( Sx_ );
	Sx_L_ = cholesky.matrixL();
	isValid_Sx_L_ = true;
      }
    
      int n = Sx_L_.cols();
      for(int i = 0; i < n; i++){
	indep_noise(i) = genGaussian_();
      }
      x_ += Sx_L_ * indep_noise;

    }

  protected:

    Vec x_; /**< \brief State */
    TimeStamp t_; /**< \brief timestamp */
    
  private:

    Mat Sx_; /**< Covariance */
    Mat Sx_inv_; /**< Inverse covariance */
    bool isValid_Sx_inv_; /**< Inverse covariance is up to date */
    double Sx_det_; /**< Determinant of Sx_ */
    double gaussian_pdf_factor_; /**< \f[ \sqrt{ (2\pi)^n)|\Sigma| } \f]*/
    bool isValid_Sx_det_; /**< Determinant of Sx_ is up to date */
    Mat Sx_L_; /**< Lower triangular part of Cholesky decomposition on Sx_ */
    bool isValid_Sx_L_; /**< Lower triangular part of Cholesky decomposition on Sx_ is up to date */

    Vec e_; /**< temporary */

    /** normal distribution random number generator */ 
    static ::boost::variate_generator< ::boost::mt19937, 
				       ::boost::normal_distribution<double> > genGaussian_;

    /** \brief Dimensionality check during initialization */
    void dimCheck(){
      assert(nDim > 0);
    }

  };

  template<unsigned int nDim>
  ::boost::variate_generator< ::boost::mt19937, 
			      ::boost::normal_distribution<double> >
  RandomVec<nDim>::genGaussian_ =
    ::boost::variate_generator< ::boost::mt19937, 
				::boost::normal_distribution<double> >
    (::boost::mt19937(rand()), ::boost::normal_distribution<double>());

} // namespace rfs

#endif
