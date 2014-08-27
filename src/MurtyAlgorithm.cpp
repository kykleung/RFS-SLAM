/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung
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

#include <cmath>
#include <stdio.h>
#include "MurtyAlgorithm.hpp"

namespace rfs
{

/// Node class implementation 

Node::Node(int n_children_exp):parent_(NULL), nChildren_(0){
  if(n_children_exp > 0)
    children_.resize(n_children_exp);
}
  
Node::~Node(){
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

Node* Node::getParent(){
  return parent_;
}

void Node::getChildren(std::vector<Node*> &children){
  children.clear();
  for( std::list<Node*>::iterator it = children_.begin(); it != children_.end(); it++ )
    children.push_back(*it);
}

int Node::getChildrenCount(){
  return nChildren_;
}

void Node::addChild(Node* child){
  child->parent_ = this;
  children_.push_front(child);
  nChildren_++;
}

/// MurtyNode class implementation

MurtyNode::MurtyNode(int id, int n_children_exp) : Node(n_children_exp),
						   id_(id),
						   assignment_(NULL), 
						   score_(0){}

  
MurtyNode::~MurtyNode(){
  if( assignment_ != NULL )
    delete[] assignment_;
  assignment_ = NULL;
}

bool MurtyNode::operator< (const MurtyNode& other){
  if( score_ < other.score_ )
    return true;
  return false;
}

void MurtyNode::setAssignment(int* a, double s){
  if(assignment_ != NULL)
    delete[] assignment_;
  assignment_ = a;
  score_ = s;
}

double MurtyNode::getScore(){
  return score_;
}

int* MurtyNode::getAssignment(){
  return assignment_;
}

int MurtyNode::getPartitionID(){
  return id_;
}

Murty::Murty(double** C, unsigned int n, double bigNum) : k_(0), n_(n), bigNumber_(bigNum){
  C_ = C;
  C_t_ = new double* [n_];
  for( int i = 0; i < n_; i++){
    C_t_[i] = new double [n_];
  }
  root_ = new MurtyNode(0, n_);
  bestScore_ = 0;
  ideal_nC_ = n;
  ideal_nR_ = n;
}
  
Murty::~Murty(){
  for(int i = 0; i < n_; i++)
    delete[] C_t_[i];
  delete[] C_t_;
  delete root_;
}

void Murty::setIdealBlock(unsigned int nR, unsigned int nC){
  ideal_nC_ = nC;
  ideal_nR_ = nR;
  if(ideal_nC_ > n_)
    ideal_nC_ = n_;
  if(ideal_nR_ > n_)
    ideal_nR_ = n_;
}

double Murty::getBestScore(){
  return bestScore_;
}

int Murty::findNextBest( int* &assignment, double* score){
    
  if( k_ == 0 ){
    int* a = new int[n_];
    double s; 
    lam_.run(C_, n_, a, &s);
    root_->setAssignment(a, s);
    k_++;
    pq.push(root_);
    assignment = a;
    *score = s;
    bestScore_ = s;
    return k_;
  }

  // pop the highest-score node
  if(pq.empty()){
    //printf("No more solutions available\n");
    assignment = NULL;
    *score = 0;
    return -1;
  }
  MurtyNode* parent = pq.top();
  int parent_partition = parent->getPartitionID();
  pq.pop();

  int* a_parent = parent->getAssignment();
  bool Z_used_by_real_m[n_];
  for(int i = 0; i < n_; i++){
    Z_used_by_real_m[i] = false;
  }
  for(int i = 0; i < ideal_nR_; i++){
    Z_used_by_real_m[ a_parent[i] ] = true;
  }

  // partition this node (into n_ - 1 - parent_partition parts) -- make children nodes
  //printf("\n************************\n");
  for(int n = parent_partition; n < n_ - 1; n++ ){ 

    if(n >= ideal_nR_ && a_parent[n] >= ideal_nC_ && Z_used_by_real_m[n - ideal_nR_]){
      // printf("Skipping partition due to redundant assignment");
      continue;
    }
      
    //printf("\nPartition %d\n", n);
    MurtyNode* p = new MurtyNode(n, n_);
    parent->addChild(p);
    int* a = new int[n_];
    bool freeCol[n_];
    for(int i = 0; i < n_; i++){
      freeCol[i] = true;
    }
    double assignmentFixedScore = 0;

    // Copy all the fixed assignments from the parent node
    //printf("Must have: ");
    for(int i = 0; i < parent_partition; i++){
      a[i] = a_parent[i];
      freeCol[a[i]] = false;
      assignmentFixedScore += C_[i][a[i]];
      //printf("(%d,%d)[%f] ", i, a[i], C_[i][a[i]]);
    }

    // assignments that have to be included in this partition
    for(int i = parent_partition; i < n; i++){
      a[i] = a_parent[i];
      freeCol[a[i]] = false;
      assignmentFixedScore += C_[i][a[i]];
      //printf("(%d,%d)[%f] ", i, a[i], C_[i][a[i]]);
    } 
    //printf("\n");

    // build up the cost matrix to solve, which may be smaller than the original cost matrix due to fixed assignments
    int nFreeAssignments = n_ - n;
    int rowRemap[nFreeAssignments]; // find entry in original matrix from reduced matrix
    int rowRemapR[n_]; // find extry in reduced matrix from original matrix
    int colRemap[nFreeAssignments];
    int colRemapR[n_];
    int nFreeCols = 0;

    for(int i = 0; i < nFreeAssignments; i++){
      rowRemap[i] = n + i;
      rowRemapR[n+i] = i;
      //printf("Remap row %d to %d\n", rowRemap[i], i);
    }
    for(int j = 0; j < n_; j++){
      if(freeCol[j]){
	colRemap[nFreeCols] = j;
	colRemapR[j] = nFreeCols;
	//printf("Remap col %d to %d\n", colRemap[nFreeCols], nFreeCols);
	nFreeCols++; 
      }
    } // at this point nFreeCol should equal nFreeAssignments
    for(int i = 0; i < nFreeAssignments; i++ ){
	
      for(int j = 0; j < nFreeAssignments; j++ ){
	//printf("Remap [%d,%d] to [%d,%d]\n", original_i, original_j, i, j);
	C_t_[i][j] = C_[rowRemap[i]][colRemap[j]];
      }
    }

    // Set all the assignments that cannot occur
    //printf("Must not have: ");
    MurtyNode* current = p;
    MurtyNode* next;
    do{
      int currentPart = current->getPartitionID();
      if( currentPart < p->getPartitionID() ){
	break;
      }
      next = static_cast< MurtyNode* > ( current->getParent() );
      int* nextAssignment = next->getAssignment();
      int doNotAssign_i = rowRemapR[currentPart];
      int doNotAssign_j = colRemapR[nextAssignment[currentPart]];
      C_t_[doNotAssign_i][doNotAssign_j] = -bigNumber_;	
      //printf("(%d,%d)=>(%d,%d) ", currentPart, nextAssignment[currentPart],doNotAssign_i,doNotAssign_j);
      current = next;
	
    }while (current != root_ );
    //printf("\n");
    bool solutionPossible = false;
    int constraintRow = rowRemapR[p->getPartitionID()];
    for(int j = 0; j < nFreeAssignments; j++){
      if( C_t_[constraintRow][j] != -bigNumber_){
	solutionPossible = true;
	break;
      }
    }

    // Solve best-linear assignment for the n-th partition
    if(solutionPossible){
      int aTmp[nFreeAssignments];
      double s = 0;
      /*printf("Solving best linear assignment for cost matrix:\n");
	for(int i = 0; i < nFreeAssignments; i++ ){
	for(int j = 0; j < nFreeAssignments; j++ ){
	printf("%f   ", C_t_[i][j]);
	}
	printf("\n");
	}*/
      lam_.run(C_t_, nFreeAssignments, aTmp, &s);
      
      double s_more_accurate = 0;
      for(int i = 0; i < nFreeAssignments; i++){
	int i_actual = rowRemap[i];
	int j_actual = colRemap[aTmp[i]];
	a[i_actual] = j_actual;
	s_more_accurate += C_[i_actual][a[i_actual]];
      }
      /*printf("Assignment:\n");
	for(int i = 0; i < n_; i++){
	printf("x[%d] ----- y[%d]\n", i, a[i]);
	}*/
      s_more_accurate += assignmentFixedScore;
      // printf("Score: %f\n", s);
      
      // Insert the solution in the priority queue and set it as a child of the highest-score node
      p->setAssignment(a, s_more_accurate);
      pq.push(p);
    }else{
      //printf("No solution possible with this partition\n");
      for(int i = 0; i < n_; i++){
	a[i] = -1;
      }
      p->setAssignment(a, -1);
    }
      
  } // end for all partitions

    // All new partitions have been inserted, now look for the highest score
  if(pq.empty()){
    //printf("No more solutions available\n");
    assignment = NULL;
    *score = 0;
    return -1;
  }
  MurtyNode* node_hi = pq.top();
  assignment = node_hi->getAssignment();
  *score = node_hi->getScore();
    
  k_++;

  return k_;
}

}
