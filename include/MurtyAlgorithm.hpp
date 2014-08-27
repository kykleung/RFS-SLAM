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

#ifndef MURTY_ALGORITHM_HPP
#define MURTY_ALGORITHM_HPP

#include <cstddef>
#include <list>
#include <queue> 
#include "HungarianMethod.hpp"

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
  Node(int n_children_exp = -1);
  
  /** Destructor, virtual so that the children of derived objects gets destroyed properly */
  virtual ~Node();

  /** 
   * Get a pointer to the parent node
   * \return pointer to parent, or NULL if this node has no parent
   */
  Node* getParent();

  /**
   * Get pointers to the children
   * \param[out] children a vector containing pointers to children nodes
   * \return a pointers to a vector
   */ 
  void getChildren(std::vector<Node*> &children);

  /**
   * Get the number of children that this node has
   * \return number of children
   */
  int getChildrenCount();

  /**
   * Add a child to this node. 
   * \note The parent node will take care of memory deallocation of children nodes
   * \param[in] child pointer to child node
   */
  void addChild(Node* child);

private:

  std::list<Node*> children_;
  int nChildren_;
  Node* parent_;
};



/**
 * \class MurtyNode
 * A node used in the Murty class for k-best linear assignment
 * \brief A node used in the Murty class for k-best linear assignment
 */
class MurtyNode : public Node
{
public:

  /** Constructor 
   * \param[in] id The partition id number for this node. Default is -1
   * \param[in] n_children_exp the number of children this node is expected to have,
   * used for more efficient memory allocation. Not mandatory.
   */
  MurtyNode(int id = -1, int n_children_exp = -1);

  /** Destructor */
  ~MurtyNode();

  /** Comparison operator for sorting */
  bool operator< (const MurtyNode& other);

  void setAssignment(int* a, double s);

  double getScore();

  int* getAssignment();

  int getPartitionID();

private:

  int id_; 
  int* assignment_;
  double score_;

};


/**
 * \class MurtyNodeCompare
 * A comparison class used for the priorty queue in Murty
 * \brief A comparison class used for the priorty queue in Murty
 */
class MurtyNodeCompare{
public:
  bool operator()(MurtyNode* p1, MurtyNode* p2){
    if( p1->getScore() < p2->getScore() )
      return true;
    return false;
  }
};


/**
 * \class Murty
 * Murty's k-best linear assignment algorithm
 * Given the partitioning of a linear assignment problem, and the
 * current best possible assignment,  Murty's algorithm will find 
 * the next best assignment
 * \brief Murty's k-best linear assignment algorithm
 */ 
class Murty{

public:
   
  /** constructor 
   *  \param[in] C square score matrix from which we will find the best assignments.
   *  \param[in] n dimension of C
   *  \param[in] bigNum A number whose negative is used to replace entries in the cost matrix that cannot be used for assignment
   */ 
  Murty(double** C, unsigned int n, double bigNum = 10000);
  
  /** destructor */
  ~Murty();

  /**
   * Get the best score.
   * \note The best score = 0 when object is first initiated.
   * \return the best score
   */
  double getBestScore();

  /** Find the next best assignment 
   *  \param[out] assignment pointer to array of assignment (no need to allocate/deallocate )
   *  \param[out] score the assignment score 
   *  \return k, the k-best solution index, or -1 if no more solutions are available
   */
  int findNextBest( int* &assignment, double* score);

  /** This is useful where we need to consider cases in whcih it is possible for elements
   *  to be unassigned. For example, in pairing landmarks with measurements, landmarks may be 
   *  misdetected, and measurements may be false alarm. To use this function, the origin (ideal)
   *  m x n assignment matrix should be augmented to be a (m+n) x (m+n) matrix, such that the
   *  ideal matrix is in the upper left. The upper right augmented block should represent miss-detections 
   *  and the lower left augmented block should represent false alarms. 
   *  \param[in] nR number of rows in the ideal block
   *  \param[in] nC number of columns in the ideal block
   */
  void setIdealBlock(unsigned int nR, unsigned int nC);

private:
  
  MurtyNode* root_; /**< The best linear assignment */
  int k_; /**< number of solutions found */
  double** C_; /**< the nxn score matrix */
  double** C_t_; /**< temporary score matrix */
  unsigned int n_; /**< dimension of C_ and C_t_*/
  HungarianMethod lam_; /**< Linear assignment method */
  std::priority_queue<MurtyNode*, std::vector<MurtyNode*>, MurtyNodeCompare> pq; /**< Priority queue for easy retrieval of node with highest score */
  double bestScore_;
  double bigNumber_;
  unsigned int ideal_nC_; /**< dimension of ideal assignment block */
  unsigned int ideal_nR_; /**< dimention of ideal assignment block */

};

}

#endif

