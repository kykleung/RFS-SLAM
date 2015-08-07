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

  //------ Node Implementation -----//

  Node::Node():parent_(NULL){}
  
  Node::~Node(){}

  Node* Node::getParent() const{
    return parent_;
  }

  Node* Node::getChildren(uint idx) const{
    if(idx >= getChildrenCount())
      return NULL;
    return children_[idx].get();
  }

  uint Node::getChildrenCount() const{
    return children_.size();
  }

  void Node::addChild(Node* child){
    child->parent_ = this;
    children_.push_back( Node::NodePtr(child) );
  }

  void Node::addChild(Node &child){
    child.parent_ = this;
    children_.push_back( Node::NodePtr( new Node(child) ) );
  }


  //----- MurtyNode Implementation -----//

  MurtyNode::MurtyNode(int id, int n_children_exp) : id_(id), 
						     s_(0){

    if(n_children_exp > 0){
      children_.reserve( n_children_exp );
    }
  }

  
  MurtyNode::~MurtyNode(){}

  bool MurtyNode::operator< (const MurtyNode& other){
    if( s_ < other.s_ )
      return true;
    return false;
  }

  void MurtyNode::setAssignment(MurtyNode::Assignment &a, double s){
    a_ = a;
    s_ = s;
  }

  double MurtyNode::getScore() const{
    return s_;
  }

  MurtyNode::Assignment MurtyNode::getAssignment() const{
    return a_;
  }

  int MurtyNode::getPartitionID() const{
    return id_;
  }



  Murty::Murty(double** C, unsigned int n, double bigNum) : k_(0), 
							    n_(n), 
							    bigNumber_(bigNum),
							    root_( new MurtyNode(0) ),
							    realAssign_nR_(n),
							    realAssign_nC_(n)
  {
    C_ = C;
    C_t_ = new double* [n_];
    for( int i = 0; i < n_; i++){
      C_t_[i] = new double [n_];
    }
    bestScore_ = 0;
  }
  
  Murty::~Murty(){
    for(int i = 0; i < n_; i++)
      delete[] C_t_[i];
    delete[] C_t_;
  }

  void Murty::setRealAssignmentBlock(unsigned int nR, unsigned int nC){
    realAssign_nC_ = nC;
    realAssign_nR_ = nR;
    if(realAssign_nC_ > n_)
      realAssign_nC_ = n_;
    if(realAssign_nR_ > n_)
      realAssign_nR_ = n_;
  }

  double Murty::getBestScore() const{
    return bestScore_;
  }

  int Murty::findNextBest( Assignment &assignment, double &score){

    // For notation, assume Z are the elements represented by the rows of the cost matrix
    // and M are the elements represented by the columns of the cost matrix
  
    //----- Calling for the first time -----//
    if( k_ == 0 ){
      Assignment a(new int[n_]);
      double s; 
      lam_.run(C_, n_, a.get(), &s);
      root_->setAssignment(a, s);
      k_++;
      pq.push(root_.get());
      assignment = a;
      score = s;
      bestScore_ = s;
      return k_;
    }

    //----- Pop the highest-score node -----//
    if(pq.empty()){ // no more solution available
      assignment.reset();
      score = 0;
      return -1;
    }
    MurtyNode* parent = pq.top();
    int parent_partition = parent->getPartitionID();
    pq.pop();

    Assignment a_parent = parent->getAssignment();
    bool Z_pairedToReal_M[n_];
    for(int i = 0; i < n_; i++){
      Z_pairedToReal_M[i] = false;
    }
    for(int i = 0; i < realAssign_nR_; i++){
      Z_pairedToReal_M[ a_parent[i] ] = true;
    }

    //----- partition this node (into  - 1 - parent_partition_id parts) -----//
    // printf("\n************************\n");
    int partitionMax = realAssign_nR_;
    int nNewPartitions = realAssign_nR_ - parent_partition + 1;
    if( realAssign_nR_ == n_ ){
      partitionMax = n_ - 1;
      nNewPartitions = n_ - parent_partition - 1;
    }
    // printf("Parent partition: %d, Make %d new partitions\n", parent_partition, nNewPartitions);
    for(int n = parent_partition; n < partitionMax; n++ ){ 
      
      // printf("\nPartition %d\n", n);
      MurtyNode* p = new MurtyNode(n, n_);
      parent->addChild(p);
      Assignment a(new int[n_]);
      bool freeCol[n_];
      for(int i = 0; i < n_; i++){
	freeCol[i] = true;
      }
      double assignmentFixedScore = 0;

      //----- Identify fixed assignments -----//
      // printf("Must have: ");
      // Copy all the fixed assignments from the parent node
      for(int i = 0; i < parent_partition; i++){
	a[i] = a_parent[i];
	freeCol[a[i]] = false; // M[i] has been paired with Z[ a[i] ]
	assignmentFixedScore += C_[i][a[i]];
	// printf("(%d,%d)[%f] ", i, a[i], C_[i][a[i]]);
      }
      // assignments that have to be included in this partition
      for(int i = parent_partition; i < n; i++){
	a[i] = a_parent[i];
	freeCol[a[i]] = false;
	assignmentFixedScore += C_[i][a[i]];
	//printf("(%d,%d)[%f] ", i, a[i], C_[i][a[i]]);
      } 
      // printf("\n");

      //----- build up the cost matrix to solve from non-fixed assignemnts-----// 
      int nFreeAssignments = n_ - n;
      int rowRemap[nFreeAssignments]; // find entry in original matrix from reduced matrix
      int rowRemapR[n_];              // find entry in reduced matrix from original matrix (reverse remap)
      int colRemap[nFreeAssignments];
      int colRemapR[n_];
      int nFreeCols = 0;

      for(int i = 0; i < nFreeAssignments; i++){
	rowRemap[i] = n + i;
	rowRemapR[n+i] = i;
	// printf("Remap row %d to %d\n", rowRemap[i], i);
      }
      for(int j = 0; j < n_; j++){
	if(freeCol[j]){
	  colRemap[nFreeCols] = j;
	  colRemapR[j] = nFreeCols;
	  // printf("Remap col %d to %d\n", colRemap[nFreeCols], nFreeCols);
	  nFreeCols++; 
	}
      } // at this point nFreeCol should equal nFreeAssignments
      for(int i = 0; i < nFreeAssignments; i++ ){
	for(int j = 0; j < nFreeAssignments; j++ ){
	  C_t_[i][j] = C_[rowRemap[i]][colRemap[j]];
	}
      }

      //----- Set all the assignments that cannot occur (negative constraints) ----->
      // printf("Must not have: ");
      MurtyNode* current = p;
      MurtyNode* next;
      do{
	int currentPart = current->getPartitionID();
	next = static_cast< MurtyNode* > ( current->getParent() );
	Assignment nextAssignment = next->getAssignment();
	int doNotAssign_i = rowRemapR[currentPart];
	int doNotAssign_j = colRemapR[nextAssignment[currentPart]];
	C_t_[doNotAssign_i][doNotAssign_j] = -bigNumber_;	
	if( doNotAssign_j >= realAssign_nC_ ){ // check if doNotAssign_j is a real assignment to avoid duplicate assignments
	  for(int y = 0; y < nFreeAssignments; y++ ){
	    if( y >= realAssign_nC_ ){
	      C_t_[doNotAssign_i][y] = -bigNumber_;
	    }
	  }
	}
	// printf("(%d,%d)=>(%d,%d) ", currentPart, nextAssignment[currentPart],doNotAssign_i,doNotAssign_j);
	current = next;	
      }while (current != root_.get() && current->getPartitionID() >= p->getPartitionID());
      // printf("\n");
    
      //----- Will we ever end up with a row that can't be assigned to any column? -----//
      //----- This may not be necessary ------///
      bool solutionPossible = false;
      int constraintRow = rowRemapR[p->getPartitionID()];
      for(int j = 0; j < nFreeAssignments; j++){
	if( C_t_[constraintRow][j] != -bigNumber_){
	  solutionPossible = true;
	  break;
	}
      }

      //----- Solve best-linear assignment for the partition -----//
      if(solutionPossible){
	Assignment aTmp(new int[nFreeAssignments]);
	double s = 0;

	/*printf("Solving best linear assignment for cost matrix:\n");
	for(int i = 0; i < nFreeAssignments; i++ ){
	for(int j = 0; j < nFreeAssignments; j++ ){
	printf("%f   ", C_t_[i][j]);
	}
	printf("\n");
	}*/

	lam_.run(C_t_, nFreeAssignments, aTmp.get(), &s);
      
	// Hungarian method performs operations directly on cost matrix
	// and may have introduced small numerical differences to score.
	// Therefore we calculate score from the original cost matrix.
	double s_more_accurate = 0; 
	for(int i = 0; i < nFreeAssignments; i++){
	  int i_actual = rowRemap[i];
	  int j_actual = colRemap[aTmp[i]];
	  a[i_actual] = j_actual;
	  s_more_accurate += C_[i_actual][a[i_actual]];
	}
	s_more_accurate += assignmentFixedScore;

	/*printf("Assignment:\n");
	for(int i = 0; i < n_; i++){
	  printf("x[%d] ----- y[%d]\n", i, a[i]);
	}
	printf("Score: %f\n", s_more_accurate);*/
      
	// Insert the solution in the priority queue
	p->setAssignment(a, s_more_accurate);
	pq.push(p);
      }else{ // no solution possible... again, do we ever get to this case?
	//printf("No solution possible with this partition\n");
	/*for(int i = 0; i < n_; i++){
	  a[i] = -1;
	  }
	  p->setAssignment(a, -1);*/
      }
      
    } // end for all partitions

    //----- All new partitions have been inserted, now look for the highest score -----//
    if(pq.empty()){
      assignment.reset();
      score = 0;
      return -1;
    }
    MurtyNode* node_hi = pq.top();
    assignment = node_hi->getAssignment();
    score = node_hi->getScore();
    k_++;
    return k_;
  }

}
