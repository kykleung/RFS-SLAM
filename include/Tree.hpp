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

#ifndef TREE_HPP
#define TREE_HPP

#include <cstddef>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <iostream>

typedef unsigned int uint;

namespace rfs{

/** 
 * \class TreeNode
 * \brief A generic tree node, defined using curiously recurring template pattern (CRTP)
 * \tparam Derived the derived class
 */
template<class Derived>
class TreeNode : 
    public boost::enable_shared_from_this<Derived>{

public:

  typedef boost::shared_ptr<Derived> DerivedPtr;
  typedef boost::weak_ptr<Derived> DerivedWeakPtr;

  /** 
   * \brief Constructor 
   * \param[in] n_children_exp the number of children this node is expected to have,
   * used for more efficient memory allocation. Not mandatory.
   */
  TreeNode(size_t n_children_exp = 0);
  
  /** \brief Destructor, virtual becuase other classes will be derived from this class*/
  virtual ~TreeNode();

  /** 
   * \brief Get a pointer to the parent node
   * \return pointer to parent, or NULL if this node has no parent
   */
  DerivedPtr getParent();
  
  /**
   * \brief Get pointers to a child
   * \param[in] idx child index
   * \return a pointer to the child
   */ 
  DerivedPtr getChild(size_t idx);

  /**
   * \brief Get the number of children that this node has
   * \return number of children
   */
  size_t getChildrenCount();

  /**
   * \brief Create a new child for this node. 
   * \return pointer to new child node
   */
  DerivedPtr addChild();

  /**
   * \brief Replace a current child with a new one
   * \param[in] idx Child index
   * \param[in] c new child
   * \return True if successful
   */
  bool replaceChild(uint idx, const DerivedPtr& c);

  /**
   * \brief Add a new child to this node.
   * \param child pointer to child
   */
  void addChild(DerivedPtr child);

  /**
   * \brief Remove all children
   */
  void removeChildren();

protected:

  std::vector< boost::shared_ptr<Derived> > children_; /**< \brief children of this node */
  DerivedWeakPtr parent_; /**< \brief parent of this node */
};

}

// Implementation

namespace rfs{

  template<class Derived>
  TreeNode<Derived>::TreeNode(size_t n_children_exp){
    if(n_children_exp > 0)
      children_.reserve(n_children_exp);
  }
  
  template<class Derived>
  TreeNode<Derived>::~TreeNode(){
    children_.clear();
  }
  
  template<class Derived>
  typename TreeNode<Derived>::DerivedPtr TreeNode<Derived>::getParent(){
    return parent_.lock();
  }
  
  template<class Derived>
  typename TreeNode<Derived>::DerivedPtr TreeNode<Derived>::getChild(size_t idx){
    if(idx < children_.size())
      return children_[idx];
    else
      return TreeNode<Derived>::DerivedPtr();
  }

  template<class Derived>
  size_t TreeNode<Derived>::getChildrenCount(){
    return children_.size();
  }

  template<class Derived>
  typename TreeNode<Derived>::DerivedPtr TreeNode<Derived>::addChild(){
    
    DerivedPtr child(new Derived);
    child->parent_ = this->shared_from_this();
    children_.push_back(child);
    return child;
  }

  template<class Derived>
  void TreeNode<Derived>::addChild(boost::shared_ptr<Derived> child){
    child->parent_ = this->shared_from_this();
    children_.push_back(child);
    return;
  }

  template<class Derived>
  bool TreeNode<Derived>::replaceChild(uint idx, const TreeNode<Derived>::DerivedPtr& c){
    
    if( idx >= this->getChildrenCount() )
      return false;
    this->children_[idx]->parent_.reset();
    this->children_[idx] = c;
    this->children_[idx]->parent_ = this->shared_from_this();
    return true;
  }

  template<class Derived>
  void TreeNode<Derived>::removeChildren(){
    for(int i = 0; i < children_.size(); i++){
      children_[i]->parent_ = DerivedPtr();
    }
    children_.clear();
  }

}



#endif
