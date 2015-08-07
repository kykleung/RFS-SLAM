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

#ifndef SPATIAL_INDEX_TREE_HPP
#define SPATIAL_INDEX_TREE_HPP

#include "SpatialIndexBox.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace rfs{

  /**
   * \class SpatialIndexTree
   * \tparam nDim number of dimensions
   * \tparam DataType type of object stored / associated with a node. This should be derived from the RandomVec class.
   * \brief Spatial indexing data structure (Octree when nDim = 3, Quadtree when nDim = 2)
   * \author Keith Leung
   */
  template <unsigned int nDim = 3, class DataType = RandomVec<nDim> >
  class SpatialIndexTree
  {
  public:

    typedef Box<nDim, DataType> TreeBox;
    typedef boost::shared_ptr<TreeBox> TreeBoxPtr;
    typedef Eigen::Matrix<double, nDim, 1> Pos;
    typedef boost::shared_ptr<DataType> DataPtr;

    /** 
     * \brief Constructor
     * \param minBoxSize minimal geometric size of a box 
     */
    SpatialIndexTree(double minBoxSize = 1);
    
    /** \brief Destructor */
    ~SpatialIndexTree();
    
    /**
     * \brief Add data to the tree. 
     *
     * This will use get() from DataType's base RandomVec class to obtain the position
     * from the first nDim elements.
     * \param[in] data Pointer to the data
     * \return box to which the data was added to
     */
    TreeBoxPtr addData(const DataPtr &data);

    /**
     * \brief Remove data from the tree
     * \param[in] data Pointer to data to be removed
     * \return True if successful
     */
    bool removeData(DataPtr data);

    /**
     * \brief Get the total number of data points stored in the tree
     * \return number of data points
     */
    uint getDataSize();

    /**
     * \brief Get all the data points within the query box
     * \param[in] qbox_min the min corner of the query box
     * \param[in] qbox_max the max corner of the query box
     * \param[out] query_result The result vector
     */
    void queryDataInBox(const Pos &qbox_min, const Pos &qbox_max, std::vector<DataPtr> &query_result);

    /**
     * \brief Find the closest data point to query point
     * \param[in] qp query point
     * \return Pointer to closest data point
     */
    DataPtr queryClosestPt(const Pos &qp);

    /**
     * \brief Clear the tree
     */
    void clear();

    /**
     * \brief Export tree structure to ascii file for visualization
     * \param[in] filename Filename
     */
    void exportASCII(std::string filename);
    

  protected:

    /**
     * \brief Divide a box into smaller boxes
     * \param[in] b Pointer to box
     * \return True if successful
     */
    bool branch(TreeBoxPtr b);

    /**
     * \brief Get the position association with a data point
     * \param[in] data Data point
     * \return position
     */
    Pos getDataPos(DataPtr d);

    /**
     * \brief Search for the smallest box that contains a point
     * \param[in] p query point
     * \param[in] b Box to start searching from. Default is the root of the tree (i.e., biggest bounding box)
     * \return pointer to box. NULL pointer if box cannot be found.
     */
    TreeBoxPtr search(const Pos &p, TreeBoxPtr b = TreeBoxPtr() );

    TreeBoxPtr root_; /**< \brief root of the tree */
    double minBoxSize_; /**< \brief smallest box size allowed */

  };

  /* Implementation */

  template<unsigned int nDim, class DataType>
  SpatialIndexTree<nDim, DataType>::SpatialIndexTree(double minBoxSize) : minBoxSize_(minBoxSize) {

    root_ = TreeBoxPtr(new TreeBox( minBoxSize*2, Pos::Ones() * -minBoxSize));
    branch(root_);
  }

  template<unsigned int nDim, class DataType>
  SpatialIndexTree<nDim, DataType>::~SpatialIndexTree(){}

  template<unsigned int nDim, class DataType>
  typename SpatialIndexTree<nDim, DataType>::TreeBoxPtr SpatialIndexTree<nDim, DataType>::addData(const DataPtr &data){
    
    // Data position
    Pos x = getDataPos(data);
    TreeBoxPtr b = search(x);

    // Check to see if we need to replace the current root with a new one with bigger bounding box
    while(b.get() == NULL){
      
      TreeBoxPtr root_old = root_;
      
      root_.reset( new TreeBox( root_old->size_*2, Pos::Ones() * -root_old->size_) );
      branch(root_); 

      // Children of root_old need to become children of the new branches of root_
      for(int i = 0; i < root_old->getChildrenCount(); i++){
	for(int c = 0; c < root_->getChildrenCount(); c++){

	  if(root_->getChild(c)->isInside( root_old->getChild(i)->getPos(TreeBox::POS_CENTER) )){
	    if(root_old->getChild(i)->getDataSize() > 0){

	      branch(root_->getChild(c));
	      for(int j = 0; j < root_->getChild(c)->getChildrenCount(); j++){
		if(root_->getChild(c)->getChild(j)->isInside( root_old->getChild(i)->getPos(TreeBox::POS_CENTER)) ){
		  root_->getChild(c)->replaceChild(j, root_old->getChild(i));
		  break;
		}
	      }
	    }

	    break;
	  } 
	 
	}
      }

      // check if data point can be placed in the new root, now with a bigger bounding box;
      b = search(x);
    }
    
    if(b->size_ == minBoxSize_ || b->getDataSize() == 0 ){ // b is the smallest box allowed, or it is empty
     
      b->addData(data);

    }else{ // check if adding data point will divide b

      while(b->getDataSize() > 0){

	Pos y = getDataPos( b->getData(0) ); // Get position of any arbitrary data point from b

	// Determine the smallest box that this datapoint could go in
	Pos smallest_box_min;
	for(int d = 0; d < nDim; d++){
	  smallest_box_min[d] = floor( y[d]/minBoxSize_ ) * minBoxSize_;
	}
	TreeBox smallest_box(minBoxSize_, smallest_box_min);

	// check if we need branch box b to index the new data point
	if(!smallest_box.isInside(x)){

	  if(branch(b)){

	    // place existing data points of b into children boxes
	    int n_max = b->getDataSize();
	    for(int n = 0; n < n_max; n++){
	      DataPtr d = b->getData(n);
	      TreeBoxPtr c = search( getDataPos(d), b);
	      c->addData(d);
	    }
	    // remove old copies of data in b
	    b->removeData();
	  
	    // determine which child box we should add data to
	    b = search(x, b);

	  }else{ // cannot branch b into smaller boxes
	    break;
	  }
	
	}else{ // new data point fits in the same (smallest-sized) box as an existing datapoint
	  break;
	}

      }
      
      b->addData(data);
    }

    return b;
     
  }

  template<unsigned int nDim, class DataType>
  bool SpatialIndexTree<nDim, DataType>::removeData(SpatialIndexTree<nDim, DataType>::DataPtr data){
    
    TreeBoxPtr b = search(getDataPos(data));
    if( b->removeData(data) ){
      // See if we can simplify the tree

      TreeBoxPtr p = b->getParent();
      TreeBoxPtr c;

      if(p->getParent().get() == NULL)
	return true;

      for(int i = 0; i < p->getChildrenCount(); i++){

	if(p->getChild(i)->getChildrenCount() != 0){
	  c.reset();
	  break;
	}

	if( p->getChild(i)->getDataSize() > 0 ){
	  if( c.get() == NULL )
	    c = p->getChild(i);
	  else{
	    c.reset();
	    break;
	  }
	}
      }
      while( c.get() != NULL ){ // only 1 child has data, move everything to p and remove all children
	for(int i = 0; i < c->getDataSize(); i++){
	  p->addData( c->getData(i) ); 
	}
	p->removeChildren();
	c.reset();
	p = p->getParent();
	if( p.get() != NULL && p->getParent().get() != NULL){
	  for(int i = 0; i < p->getChildrenCount(); i++){
	    if(p->getChild(i)->getChildrenCount() != 0){
	      c.reset();
	      break;
	    }
	    if( p->getChild(i)->getDataSize() > 0 && p->getChild(i)->getChildrenCount() == 0){
	      if( c.get() == NULL )
		c = p->getChild(i);
	      else{
		c.reset();
		break;
	      }
	    }
	  }
	}
      }

      return true;
    }
    return false;

  }

  template<unsigned int nDim, class DataType>
  uint SpatialIndexTree<nDim, DataType>::getDataSize(){

    return root_->getDataSize();
  }

  template<unsigned int nDim, class DataType>
  typename SpatialIndexTree<nDim, DataType>::TreeBoxPtr 
  SpatialIndexTree<nDim, DataType>::search(const SpatialIndexTree<nDim, DataType>::Pos &p,
					   SpatialIndexTree<nDim, DataType>::TreeBoxPtr b){
    
    int nChildrenPerBox = pow(2,nDim);
    
    if( b.get() == NULL ){
      b = root_;
    } 
    if( !(b->isInside(p)) ){
      return TreeBoxPtr();
    }
    
    while( b->getChildrenCount() != 0 ){
      bool childFound = false;
      for(int c = 0; c < nChildrenPerBox; c++){
	if (b->getChild(c)->isInside(p)){
	  b = b->getChild(c);
	  childFound = true;
	  break;
	}
      }
      if(!childFound){ // if program reaches this point, something is wrong
	std::cerr << "Spatial index tree search error\n";
	throw;
      }
    }
    return b;

  }
  
  template<unsigned int nDim, class DataType>
  bool SpatialIndexTree<nDim, DataType>::branch(SpatialIndexTree<nDim, DataType>::TreeBoxPtr b){
    
    if(b->size_ == minBoxSize_){
      return false;
    }

    int nChildrenPerBox = pow(2,nDim);
    double childBoxSize = b->size_ / 2.0;
    for(int i = 0; i < nChildrenPerBox; i++){
      
      Pos x = Pos::Zero();
      for(int d = 0; d < nDim; d++){
	x[d] = ((i >> d) & 1) == 1; 
      }
      x *= childBoxSize;
      x += b->getPos();

      b->addChild( TreeBoxPtr(new TreeBox(childBoxSize, x)) );

    }

    return true;
  }

  template<unsigned int nDim, class DataType>
  typename SpatialIndexTree<nDim, DataType>::Pos 
  SpatialIndexTree<nDim, DataType>::getDataPos(SpatialIndexTree<nDim, DataType>::DataPtr d){

    return d->get().template head<nDim>();
  }

  template<unsigned int nDim, class DataType>
  void SpatialIndexTree<nDim, DataType>::queryDataInBox(const Pos &qbox_min, const Pos &qbox_max, std::vector<DataPtr> &query_result){

    // Find the smallest box in the tree that encapsulates the query box
    TreeBoxPtr b = root_;
    bool smallerBoxFound = true;
    while( smallerBoxFound == true ){
      smallerBoxFound = false;
      for(int i = 0; i < b->getChildrenCount(); i++){
	TreeBoxPtr c = b->getChild(i);
	if( c->isInside( qbox_min ) && c->isInside( qbox_max ) ){
	  b = c;
	  smallerBoxFound = true;
	  break;
	} 
      }
    }

    std::vector<TreeBoxPtr> traverseStack;
    traverseStack.push_back(b);

    while(!traverseStack.empty()){

      // pop from stack
      b = traverseStack.back();
      traverseStack.pop_back();

      // if there are children push children into stack
      if(b->getChildrenCount() > 0){
	for(int i = b->getChildrenCount() - 1; i >= 0 ; i--){	 
	  if( b->getChild(i)->getDataSize() > 0){  // but only if the children are not empty
	    traverseStack.push_back( b->getChild(i) );
	  }
	}
      }else{ // if no children, do check to see if there is data inside query box

	Pos b_max = b->getPos( TreeBox::POS_UPPER );
	Pos b_min = b->getPos( TreeBox::POS_LOWER );
	// whole box b fits in query box
	if( (qbox_max - b_max).minCoeff() >= 0 && (qbox_min - b_max).maxCoeff() <= 0 &&
	    (qbox_max - b_min).minCoeff() >= 0 && (qbox_min - b_min).maxCoeff() <= 0){
	  for(int i = 0; i < b->getDataSize(); i++){
	    query_result.push_back( b->getData(i) );
	  }
	}else{ // box b partially fix in query box
	  for(int i = 0; i < b->getDataSize(); i++){
	    Pos x = getDataPos(b->getData(i));
	    if( (qbox_max - x).minCoeff() >= 0 && (qbox_min - x).maxCoeff() <= 0 ){
	      query_result.push_back( b->getData(i) );
	    }
	  }
	}

      }

    } // while(!traverStack.empty())

  }

  template<unsigned int nDim, class DataType>
  typename SpatialIndexTree<nDim, DataType>::DataPtr 
  SpatialIndexTree<nDim, DataType>::queryClosestPt(const SpatialIndexTree<nDim, DataType>::Pos &qp){

    TreeBoxPtr b = search(qp);
    if( b.get() == NULL ){
      b = root_;
    }
    
    // Get a point d from within b
    DataPtr d;
    while(d.get() == NULL){
      
      if(b->getDataSize() > 0 ){ // b is not empty
	if(b->getChildrenCount() == 0){ // b does not have children, so b must contain data
	  d = b->getData(0);
	}else{ // b has children
	  for(int i = 0; i < b->getChildrenCount(); i++){
	    if( b->getChild(i)->getDataSize() > 0){
	      b = b->getChild(i);
	      break;
	    }
	  }
	}
      }else{ // b is empty
	b = b->getParent();
      }
    }
 
    // Construct a query box using qp and d
    Pos dp = getDataPos(d);
    double dist = (dp - qp).norm();
    Pos qbox_min = qp - Pos::Ones() * dist;
    Pos qbox_max = qp + Pos::Ones() * dist;

    // Find the smallest box that contains the query box
    b = root_;
    bool smallerBoxFound = true;
    while( smallerBoxFound == true ){
      smallerBoxFound = false;
      for(int i = 0; i < b->getChildrenCount(); i++){
	TreeBoxPtr c = b->getChild(i);
	if( c->isInside( qbox_max, qbox_min ) ){
	  b = c;
	  smallerBoxFound = true;
	  break;
	} 
      }
    }

    // from b, traverse tree to look for a closer point
    std::vector<TreeBoxPtr> traverseStack;
    traverseStack.push_back(b);

    while(!traverseStack.empty()){

      // pop from stack
      b = traverseStack.back();
      traverseStack.pop_back();

      // if there are children push children into stack if they are within the shortest found distance between qp and dp
      if(b->getChildrenCount() > 0){
	for(int i = b->getChildrenCount() - 1; i >= 0 ; i--){
	  if( b->getChild(i)->isWithinDistance(qp, dist) )
	    traverseStack.push_back( b->getChild(i) );
	}
      }else{ // if no children, go through stored data to see if a closer point dp can be found

	for(int i = 0; i < b->getDataSize(); i++){
	  double d_temp =  (getDataPos( b->getData(i) ) - qp).norm();
	  if( d_temp < dist ){
	    d = b->getData(i);
	    dp = getDataPos(d);
	    dist = d_temp;	    
	  }
	}

      }
    }
    return d;

  }

  template<unsigned int nDim, class DataType>
  void SpatialIndexTree<nDim, DataType>::clear(){

    root_.reset( new TreeBox( minBoxSize_*2, Pos::Ones() * -minBoxSize_) );
    branch(root_);
  }

  template<unsigned int nDim, class DataType>
  void SpatialIndexTree<nDim, DataType>::exportASCII(std::string filename){

    std::ofstream dataFile;
    dataFile.open(filename.c_str());
    
    TreeBoxPtr b = root_;
    std::vector<TreeBoxPtr> traverseStack;
    traverseStack.push_back(b);

    while(!traverseStack.empty()){
      
      // pop from stack
      b = traverseStack.back();
      traverseStack.pop_back();

      // draw
      Pos b_min = b->getPos(TreeBox::POS_LOWER);
      Pos b_max = b->getPos(TreeBox::POS_UPPER);
      dataFile << std::fixed << std::setprecision(3) << "BOX ";
      for(int i = 0; i < nDim; i++){
	dataFile << std::setw(15) << b_min[i];
      }
      for(int i = 0; i < nDim; i++){
	dataFile << std::setw(15) << b_max[i];
      }
      dataFile << std::setw(15) << b->getDataSize() << std::endl;

      if(b->getChildrenCount() > 0){
	for(int i = b->getChildrenCount() - 1; i >= 0 ; i--){
	  traverseStack.push_back( b->getChild(i) );
	}
      }else{
	for(int i = 0; i < b->getDataSize(); i++){
	  Pos dp = getDataPos( b->getData(i) );
	  dataFile << std::fixed << std::setprecision(3) << "PT ";
	  for(int i = 0; i < nDim; i++){
	    dataFile << std::setw(15) << dp[i];
	  }
	  for(int i = 0; i < nDim; i++){
	    dataFile << std::setw(15) << dp[i];
	  }
	  dataFile << std::setw(15) << 0 << std::endl;
	}
      }
      
    }
    
    dataFile.close();
    
  }


}
#endif
