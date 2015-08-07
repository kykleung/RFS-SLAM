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

#ifndef SPATIAL_INDEX_BOX_HPP
#define SPATIAL_INDEX_BOX_HPP

#include <assert.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>

#include "RandomVec.hpp"
#include "Tree.hpp"

namespace rfs{

  /**
   * \class Box
   * \tparam nDim number of dimensions
   * \tparam DataType type of object stored / associated within a box 
   * \brief Box used for quadtree / octree
   */
  template <unsigned int nDim = 3, class DataType = RandomVec<3> >
  class Box : public TreeNode< Box<nDim, DataType> >
  {
  public:

    typedef Eigen::Matrix<double, nDim, 1> Pos; /**< Location of a point */
    typedef boost::shared_ptr<DataType> DataPtr; /**< Data pointer */
    typedef boost::shared_ptr< Box<nDim, DataType> > Ptr; /**< Pointer to box */

    enum POSITION{
      POS_LOWER, POS_UPPER, POS_CENTER
    };
    
    /** 
     * \brief Constructor 
     * \param[in] length size of box
     * \param[in] pos Coordinates specifying the location of the minimum corner.
     */
    Box(double length = 1, 
	Eigen::Matrix<double, nDim, 1> pos = Eigen::Matrix<double, nDim, 1>::Zero() );

    /**
     * \brief Copy Constructor
     * \param[in] other The other Box
     */
    Box(const Box &other);

    /** \brief Destructor */
    virtual ~Box();

    /** 
     * \brief Assignment operator 
     * \warn It is the responsibility of the caller to ensure that the bounds are equal
     */
    Box& operator=(const Box& other);
    
    
    /** 
     * \brief Get the position of the box 
     * \param[in] p which part of the box 
     */
    Eigen::Matrix<double, nDim, 1> getPos(POSITION p = POS_LOWER);

    /**
     * \brief Add a child box
     * \param[in] b child box to add
     */
    void addChild(const Ptr& b);

    /**
     * \brief Replace an existing child box with a new ond
     * \param[in] idx index of child to be replaced
     * \param[in] b new child box
     * \return True if successful
     */
    bool replaceChild(uint idx, const Ptr& b);

    /**
     * \brief Remove all children
     */
    void removeChildren();

    /**
     * \brief Add data to box
     * \param[in] data Data to add
     */
    void addData(const DataPtr &data);

    /**
     * \brief Get the number of data points stored in this box
     */
    unsigned int getDataSize();

    /**
     * \brief Get a data point 
     * \param[in] idx Data point index.
     * \return pointer to Data, or NULL pointer if index is invalid
     */
    DataPtr getData(unsigned int idx);

    /**
     * \brief Remove a data point according to its index, or remove all data points stored in the box
     * \info for single data point removal, the last point in the list is used to replace the erased point
     * \param[in] idx Data point index. Default of -1 will erase all data points
     * \return True if successful
     */
    bool removeData(int idx = -1);

    /**
     * \brief Remove a data point by comparing actual data stored
     * \param[in] data Data point
     * \return True if successful
     */
    bool removeData(DataPtr data);

    /**
     * \brief Check if a point is within this box
     * \param[in] p query point
     * \return True or False
     */
    bool isInside(const Pos &p);

    /**
     * \brief Check if a box is within this box
     * \param[in] b query box
     * \return True or False
     */
    bool isInside(const Ptr &b);

    /**
     * \brief Check if a box is within this box
     * \param[in] b_max query box max corner
     * \param[in] b_min query box min corner
     * \return True or False
     */
    bool isInside(const Pos &b_max, const Pos &b_min);

    /**
     * \brief Check if a point is within a distance d from the box
     * \param[in] p query point
     * \param[in] d query distance
     * \return True or False
     */
    bool isWithinDistance(Pos p, double d);

    double const size_; /**< size of the box */
   

  protected:

    Pos bound_min_; /**< \brief point defining the intersection of all minimum bounds */
    Pos bound_max_; /**< \brief point defining the intersection of all maximum bounds */

    std::vector< DataPtr > data_; /**< \brief data storage */

    unsigned int nData_; /**< \brief number of data points contained */

  };

  /* Implementation */

  template <unsigned int nDim, class DataType>
  Box<nDim, DataType>::Box(double length, Eigen::Matrix<double, nDim, 1> pos) : size_(length){

    assert(nDim > 0);

    bound_min_ = pos;
    bound_max_ = pos + Pos::Ones() * length;

    nData_ = 0;
    
  }
  
  template <unsigned int nDim, class DataType>
  Box<nDim, DataType>::Box(const Box &other) : TreeNode< Box<nDim, DataType> >(other){
    
    bound_min_ = other.bound_min_;
    bound_max_ = other.bound_max_;

    int nDataDiff = other.nData_ - nData_;
    nData_ = other.nData_;

    Ptr parent = this->getParent();
    while( parent.get() != NULL ){
      parent->nData_+= nDataDiff;
      parent = parent->getParent();
    }
  }

  template <unsigned int nDim, class DataType>
  Box<nDim, DataType>::~Box(){

    removeData(); // to update parents' nData_ count
  }

  template <unsigned int nDim, class DataType>
  Box<nDim, DataType>& Box<nDim, DataType>::operator=(const Box<nDim, DataType> &other){

    TreeNode< Box<nDim, DataType> >::operator=(other);

    bound_min_ = other.bound_min_;
    bound_max_ = other.bound_max_;

    int nDataDiff = other.nData_ - nData_;
    nData_ = other.nData_;

    Ptr parent = this->getParent();
    while( parent.get() != NULL ){
      parent->nData_+= nDataDiff;
      parent = parent->getParent();
    }
  }
 
  template <unsigned int nDim, class DataType>
  typename Box<nDim, DataType>::Pos Box<nDim, DataType>::getPos(Box<nDim, DataType>::POSITION p){
    
    switch(p){
    case POS_LOWER: return bound_min_;
    case POS_UPPER: return bound_max_;
    case POS_CENTER: return (bound_max_ - bound_min_) / 2 + bound_min_;
    }
    return bound_min_;
  }

  template <unsigned int nDim, class DataType>
  void Box<nDim, DataType>::addChild(const Ptr& b){

    TreeNode< Box<nDim, DataType> >::addChild(b);
    nData_ += b->getDataSize();
  }

  template <unsigned int nDim, class DataType>
  bool Box<nDim, DataType>::replaceChild(uint idx, const Ptr& b){

    if( idx >= this->getChildrenCount() )
      return false;

    int nDataDiff = b->getDataSize() - this->getChild(idx)->getDataSize();
    
    TreeNode< Box<nDim, DataType> >::replaceChild(idx, b);

    nData_ += nDataDiff;
    Ptr parent = this->getParent();
    while( parent.get() != NULL ){
      parent->nData_+= nDataDiff;
      parent = parent->getParent();
    }

    return true;
  }

  template <unsigned int nDim, class DataType>
  void Box<nDim, DataType>::removeChildren(){

    int nDataDiff = 0;
    for(int i = 0; i < this->getChildrenCount(); i++){
      nDataDiff -= this->getChild(i)->getDataSize();
    }
    nData_ += nDataDiff;
    Ptr parent = this->getParent();
    while( parent.get() != NULL ){
      parent->nData_+= nDataDiff;
      parent = parent->getParent();
    }
    TreeNode< Box<nDim, DataType> >::removeChildren();
  }

  template <unsigned int nDim, class DataType>
  bool Box<nDim, DataType>::isInside(const Box<nDim, DataType>::Pos &p){
    
    if( (bound_max_ - p).minCoeff() >= 0 && (bound_min_ - p).maxCoeff() <= 0 ){
      return true;
    }else{
      return false;
    }
  }

  template <unsigned int nDim, class DataType>
  bool Box<nDim, DataType>::isInside(const Ptr &b){

    return isInside(b->getPos(POS_UPPER), b->getPos(POS_LOWER));
  }

  template <unsigned int nDim, class DataType>
  bool Box<nDim, DataType>::isInside(const Pos &b_max, const Pos &b_min){

    if( isInside( b_max ) && isInside( b_min ) )
      return true;
    return false;
  }

  template <unsigned int nDim, class DataType>
  void Box<nDim, DataType>::addData( const Box<nDim, DataType>::DataPtr &data){
    data_.push_back(data);
    nData_++;
    Ptr parent = this->getParent();
    while( parent.get() != NULL ){
      parent->nData_++;
      parent = parent->getParent();
    }
  }

  template <unsigned int nDim, class DataType>
  unsigned int Box<nDim, DataType>::getDataSize(){
    return nData_;
  }

  template <unsigned int nDim, class DataType>
  typename Box<nDim, DataType>::DataPtr Box<nDim, DataType>::getData(unsigned int idx){
    
    if(data_.size() > idx){
      return data_[idx];
    }else{
      return DataPtr();
    }
  }

  template <unsigned int nDim, class DataType>
  bool Box<nDim, DataType>::removeData(int idx){

    uint nDataRemoved = 0;
    if(idx < 0){
      uint nData_allChildren = 0;
      for(int i = 0; i < this->getChildrenCount(); i++){
	nData_allChildren += this->getChild(i)->getDataSize();
      }
      nDataRemoved = this->nData_ - nData_allChildren;
      data_.clear();
    }else{
      nDataRemoved = 1;
      if(idx >= data_.size() ){
	return false;
      }
      data_[idx] = data_.back();
      data_.pop_back();
    }
    nData_ -= nDataRemoved;

    Ptr parent = this->getParent();
    while( parent.get() != NULL ){
      parent->nData_-= nDataRemoved;
      parent = parent->getParent();
    }

    return true;

  }

  template <unsigned int nDim, class DataType>
  bool Box<nDim, DataType>::removeData(DataPtr data){

    for(int i = 0; i < data_.size(); i++){
      if(data_[i] == data){
	return removeData(i);
      }
    }
    return false;

  }

  template <unsigned int nDim, class DataType>
  bool Box<nDim, DataType>::isWithinDistance(Pos p, double d){

    double check_box_size = size_ + 2*d;
    Box<nDim, DataType> check_box( check_box_size, getPos(POS_CENTER) - check_box_size/2 * Pos::Ones() );
    return(check_box.isInside(p));
	     
  }

  
}

#endif

