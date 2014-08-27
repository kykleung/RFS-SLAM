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
#include <list>

namespace rfs{

/** 
 * \class Node
 * A generic tree node
 * \brief A generic tree node
 */
class Node{

public:

  /** Constructor 
   * \param[in] n_children_exp the number of children this node is expected to have,
   * used for more efficient memory allocation. Not mandatory.
   */
  Node(int n_children_exp = -1):parent_(NULL), nChildren_(0){
    if(n_children_exp > 0)
      children_.resize(n_children_exp);
  }
  
  /** Destructor, virtual so that the children of derived objects gets destroyed properly */
  virtual ~Node(){
    Node* parent = NULL;
    Node* current = this;
    do{
      if(current->nChildren_ > 0){
	parent = current;
	current = parent->children_.front();
	parent->children_.pop_front();
	parent->nChildren_--;
      }else if(current != this){
	delete current;
	current = parent;
	parent = current->parent_;
      }
    }while(parent!=NULL || current->nChildren_ > 0);
  }

  /** 
   * Get a pointer to the parent node
   * \return pointer to parent, or NULL if this node has no parent
   */
  Node* getParent(){
    return parent_;
  }
  
  /**
   * Get pointers to the children
   * \param[out] children a vector containing pointers to children nodes
   * \return a pointers to a vector
   */ 
  void getChildren(std::vector<Node*> &children){
    children.clear();
    for( std::list<Node*>::iterator it = children_.begin(); it != children_.end(); it++ )
      children.push_back(*it);
  }

  /**
   * Get the number of children that this node has
   * \return number of children
   */
  int getChildrenCount(){
    return nChildren_;
  }

  /**
   * Add a child to this node. 
   * \note The parent node will take care of memory deallocation of children nodes
   * \param[in] child pointer to child node
   */
  void addChild(Node* child){
    child->parent_ = this;
    children_.push_front(child);
    nChildren_++;
  }

private:

  std::list<Node*> children_;
  int nChildren_;
  Node* parent_;
};

}

#endif
