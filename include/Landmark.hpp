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

// Landmark classes for defining map feature state
// Keith Leung 2013

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include <algorithm>
#include "RandomVec.hpp"

namespace rfs
{

  /** 
   * \class Landmark
   * \brief An abstract class for defining landmark state
   * \tparam dim Number of dimensions
   * \tparam nVecSpaceDim The number of dimensions of the vector space that the landmark exists in. 
   * For example, 2 for 2D space (default), 3 for 3D space. By default, the first nVecSpaceDim
   * components of VecType are assumed to represent position in the vector space. All other components
   * are assumed to represent orientation.
   * \tparam DescriptorType Landmark descriptor object type
   * \author Keith Leung
   */

  template<unsigned int dim, unsigned int nVecSpaceDim = dim, class DescriptorType = int>
  class Landmark : public RandomVec<dim>
  {
  public:

    typedef typename RandomVec<dim>::Vec Vec;
    typedef typename RandomVec<dim>::Mat Mat;
    typedef typename RandomVec<dim>::Mat Cov;

    // Required because Landmark could be placed in an std::vector
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    /** \brief the dimensionality of the physical space that the landmark exists in */
    static unsigned int const nVecSpaceDim_ = nVecSpaceDim;

    /** \brief The type of information represented by each component of Vec */
    enum VecComponentType{
      VEC_COMPONENT_POSITION,
      VEC_COMPONENT_ORIENTATION,
      VEC_COMPONENT_OTHER
    };

    /** \brief Default constructor */
    Landmark(){
      for(int i = 0; i < Vec::RowsAtCompileTime; i++){
	if(i < nVecSpaceDim){
	  vecComp_[i] = VEC_COMPONENT_POSITION;
	}else{
	  vecComp_[i] = VEC_COMPONENT_ORIENTATION;
	}
      }
    }

    /** 
     * Constructor
     */
    Landmark(Vec const &x, Mat const &Sx){
      for(int i = 0; i < Vec::RowsAtCompileTime; i++){
	if(i < nVecSpaceDim){
	  vecComp_[i] = VEC_COMPONENT_POSITION;
	}else{
	  vecComp_[i] = VEC_COMPONENT_ORIENTATION;
	}
      }
      this->set(x, Sx);
    }

    /** 
     * Constructor with descriptor
     */
    Landmark(Vec &x, Mat &Sx, DescriptorType &d){
      for(int i = 0; i < Vec::RowsAtCompileTime; i++){
	if(i < nVecSpaceDim){
	  vecComp_[i] = VEC_COMPONENT_POSITION;
	}else{
	  vecComp_[i] = VEC_COMPONENT_ORIENTATION;
	}
      }
      this->set(x, Sx);
      desc_ = d;
    }

    /** Default destructor */
    ~Landmark(){};

    /** 
     * Set the type of a component of Vec 
     * \param[in] i the component index
     * \param[in] type the component type
     */
    void setComponentType(int i, VecComponentType type){
      vecComp_[i] = type;
    }

    /**
     * Get the type of a component of Vec
     * \param[in] i  the component index
     * \return the component type
     */
    VecComponentType getComponentType(int i){
      return vecComp_[i];
    }
  

    /** Set descriptor for landmark 
     *  \param[in] d descriptor
     */
    void setDescriptor(DescriptorType &d){
      desc_ = d;
    }

    /** Get descriptor for landmark 
     *  \param[in] d descriptor
     */
    void getDescriptor(DescriptorType &d){
      d = desc_;
    }

    /** 
     * Get the size of a loose bounding box that encapsulates the covariance ellipsoid.
     * This operation is based on heuristics, where a loose upper bound on the largest
     * Eigenvalue of the covariance matrix is found. The square root is then taken and
     * then multiplied by a constant (default 3). The box is axis-aligned.
     * \param[out] center pointer to an caller allocated array for storing the center of the box
     * \param[out] size pointer to an caller allocated array for storing the size of the box
     * \param[in] multiplier constant for the box size
     */
    void getBoundingBox(double* center, double* size, double multiplier = 3.0 ){

      // Get a loose bound on the largest eigenvalue of the covariance 
      Mat Sx = this->getCov();
      double rowSum[Vec::RowsAtCompileTime];
      double colSum[Vec::RowsAtCompileTime];
      for(int i = 0; i < Vec::RowsAtCompileTime; i++){
	rowSum[i] = 0;
	colSum[i] = 0;
      }
      for(int i = 0; i < Vec::RowsAtCompileTime; i++){
	for(int j = 0; j < Vec::RowsAtCompileTime; j++){
	  if(vecComp_[i] == VEC_COMPONENT_POSITION && vecComp_[j] == VEC_COMPONENT_POSITION){
	    rowSum[i] += Sx(i,j);
	    colSum[j] += Sx(i,j);
	  }
	}
      }
      double rowMax = *std::max_element(rowSum, rowSum + Vec::RowsAtCompileTime);
      double colMax = *std::max_element(colSum, colSum + Vec::RowsAtCompileTime);
      double l = 2.0 * multiplier * sqrt(fmin(rowMax, colMax));
    
      for(int i = 0; i < Vec::RowsAtCompileTime; i++){
	size[i] = l;
	center[i] = this->get(i);
      }
    }

  private:
  
    DescriptorType desc_; /**< descriptor */
    VecComponentType vecComp_[Vec::RowsAtCompileTime]; /**< description for each vector component */

  };

  typedef Landmark< 1, 1 > Landmark1d;

  typedef Landmark< 2, 2 > Landmark2d;

  typedef Landmark< 3, 3 > Landmark3d;

}

#endif
