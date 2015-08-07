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

// class for defining pose state
// Keith Leung 2013

#ifndef POSE_HPP
#define POSE_HPP

#include "RandomVec.hpp"

namespace rfs{

  /**
   * \class Pose
   * \brief A random vector with position and orientation.
   * \note The position components are always placed before orientation ones.
   * \tparam nDim Dimension of the pose vector
   * \tparam nPosDim Dimension of the position component
   * \tparam nRotDim Dimension of the rotation component
   */
  template<unsigned int nDim, unsigned int nPosDim, unsigned int nRotDim>
  class Pose : public RandomVec< nDim >
  {
  public:
    
    /** \brief the full vector */
    typedef typename RandomVec<nDim>::Vec Vec;

    /** \brief the full covariance matrix */
    typedef typename RandomVec<nDim>::Mat Cov;

    /** \brief the full covariance matrix */
    typedef typename RandomVec<nDim>::Mat Mat;

    /** \brief Position vector */
    typedef ::Eigen::Matrix<double, nPosDim, 1>PosVec; 
    
    /** \brief Position covariance matrix */
    typedef ::Eigen::Matrix<double, nPosDim, nPosDim>PosCov;

    /** \brief Orientation vector */
    typedef ::Eigen::Matrix<double, nRotDim, 1>RotVec; 

    /** \brief Orientation covariance matrix */
    typedef ::Eigen::Matrix<double, nRotDim, nRotDim>RotCov; 

    /** \brief Constructor */
    Pose()
    {
      dimCheck();
    }

    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] Sx Covariance matrix
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, Cov const &Sx, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, Sx, t) {}
    
    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] SxVec Array with diagonal entries of Covariance matrix
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, double const * const &SxVec, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, SxVec, t) {}
    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, t) {}

    /**
     * \brief Constructor with vector and covariance matrix set to zero
     * \param[in] t Timestamp
     */
    Pose(const TimeStamp &t ) : RandomVec<Vec::RowsAtCompileTime>(t) {}

    
    /** \brief Destructor */
    ~Pose(){}

    /** \brief Get the position vector 
     *  \param[out] x position vector 
     */
    void getPos(PosVec &x) const {
      x = this->x_.template head<nPosDim>();
    }

    /** \brief Get the position vector with the time stamp 
     *  \param[out] x position vector 
     *  \param[out] t timestamp 
     */
    void getPos(PosVec &x, TimeStamp &t) const{
      x = this->x_.template head<nPosDim>();
      t = this->t_;
    }

    /** \brief Get the position vector 
     *  \return position vector 
     */
    PosVec getPos() const{
      return this->x_.template head<nPosDim>();
    }

    /** \brief Get the rotation vector 
     *  \param[out] r rotation vector
     */
    void getRot(RotVec &r) const{
      r = this->x_.template segment<nRotDim>(nPosDim);
    }

    /** \brief Get the rotation vector with the time stamp 
     *  \param[out] r rotation vector
     *  \param[out] t timestamp
     */
    void getRot(RotVec &r, TimeStamp &t) const{
      r = this->x_.template segment<nRotDim>(nPosDim);
      t = this->t_;
    }

    /** \brief Get the rotation vector
     *  \return rotation vector
     */
    RotVec getRot() const{
      return this->x_.template segment<nRotDim>(nPosDim);
    }
    
  private:
    
    /** \brief Check if the template parameters for dimensionality are valid */
    void dimCheck(){
      assert(nDim >= nPosDim);
      assert(nDim >= nRotDim);
      assert(nDim == nPosDim + nRotDim);
    }

  };

  /** 
   * \brief Partial template specialization for Pose with only rotation and without position.
   */
  template<unsigned int nDim, unsigned int nRotDim>
  class Pose<nDim, 0, nRotDim> : public RandomVec< nDim >
  {
  public:

    /** \brief the full vector */
    typedef typename RandomVec<nDim>::Vec Vec;

    /** \brief the full covariance matrix */
    typedef typename RandomVec<nDim>::Mat Cov;

    /** \brief the full covariance matrix */
    typedef typename RandomVec<nDim>::Mat Mat;

    /** \brief Orientation vector */
    typedef ::Eigen::Matrix<double, nRotDim, 1>RotVec; 

    /** \brief Orientation covariance matrix */
    typedef ::Eigen::Matrix<double, nRotDim, nRotDim>RotCov; 

    /** \brief Constructor */
    Pose()
    {
      dimCheck();
    }

    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] Sx Covariance matrix
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, Cov const &Sx, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, Sx, t) {}
    
    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] SxVec Array with diagonal entries of Covariance matrix
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, double const * const &SxVec, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, SxVec, t) {}

    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, t) {}
    
    /**
     * \brief Constructor with vector and covariance matrix set to zero
     * \param[in] t Timestamp
     */
    Pose(const TimeStamp &t ) : RandomVec<Vec::RowsAtCompileTime>(t) {}
    
    /** \brief Destructor */
    ~Pose(){}

    /** \brief Get the rotation vector 
     *  \param[out] r rotation vector
     */
    void getRot(RotVec &r) const{
      r = this->x_.template head<nRotDim>();
    }

    /** \brief Get the rotation vector with the time stamp 
     *  \param[out] r rotation vector
     *  \param[out] t timestamp
     */
    void getRot(RotVec &r, TimeStamp &t) const{
      r = this->x_.template head<nRotDim>();
      t = this->t_;
    }

    /** \brief Get the rotation vector
     *  \return rotation vector
     */
    RotVec getRot() const{
      return this->x_.template head<nRotDim>();
    }
    
  private:
    
    /** \brief Check if the template parameters for dimensionality are valid */
    void dimCheck(){
      assert(nDim <= nRotDim);
    }

  };

  /** 
   * \brief Partial template specialization for Pose with only position and without rotation.
   */
  template<unsigned int nDim, unsigned int nPosDim>
  class Pose<nDim, nPosDim, 0> : public RandomVec< nDim >
  {
  public:

    /** \brief the full vector */
    typedef typename RandomVec<nDim>::Vec Vec;

    /** \brief the full covariance matrix */
    typedef typename RandomVec<nDim>::Mat Cov;

    /** \brief the full covariance matrix */
    typedef typename RandomVec<nDim>::Mat Mat;
    
    /** \brief Position vector */
    typedef ::Eigen::Matrix<double, nPosDim, 1>PosVec; 
    
    /** \brief Position covariance matrix */
    typedef ::Eigen::Matrix<double, nPosDim, nPosDim>PosCov;

    /** \brief Constructor */
    Pose(){
      dimCheck();
    }

    
    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] Sx Covariance matrix
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, Cov const &Sx, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, Sx, t) {}
    
    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] SxVec Array with diagonal entries of Covariance matrix
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, double const * const &SxVec, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, SxVec, t) {}

    /**
     * \brief Constructor 
     * \param[in] x Vector
     * \param[in] t Timestamp
     */
    Pose(Vec const &x, const TimeStamp &t = TimeStamp() ) : RandomVec<Vec::RowsAtCompileTime>(x, t) {}
    
    /**
     * \brief Constructor with vector and covariance matrix set to zero
     * \param[in] t Timestamp
     */
    Pose(const TimeStamp &t ) : RandomVec<Vec::RowsAtCompileTime>(t) {}

    /** \brief Destructor */
    ~Pose(){}

    /** \brief Get the position vector 
     *  \param[out] x position vector 
     */
    void getPos(PosVec &x) const{
      x = this->x_.template head<nPosDim>();
    }

    /** \brief Get the position vector with the time stamp 
     *  \param[out] x position vector 
     *  \param[out] t timestamp 
     */
    void getPos(PosVec &x, TimeStamp &t) const{
      x = this->x_.template head<nPosDim>();
      t = this->t_;
    }

    /** \brief Get the position vector 
     *  \return position vector 
     */
    PosVec getPos() const{
      return this->x_.template head<nPosDim>();
    }
    
  private:
    
    /** \brief Check if the template parameters for dimensionality are valid */
    void dimCheck(){
      assert(nDim < nPosDim);
    }

  };

  typedef Pose<1, 1, 0> Pose1d;
  typedef Pose<3, 2, 1> Pose2d;
  typedef Pose<6, 3, 3> Pose3d;

  typedef Pose<1, 1, 0> Position1d;
  typedef Pose<2, 2, 0> Position2d;
  typedef Pose<3, 3, 0> Position3d;


}

#endif
