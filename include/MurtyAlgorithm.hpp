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

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <cstddef>
#include <vector>
#include <queue> 
#include "HungarianMethod.hpp"

typedef unsigned int uint;

namespace rfs{

  /** 
   * \class Node
   * \brief A generic tree node.
   * \author Keith Leung
   */
  class Node{

  public:

    /** Boost scoped pointer to a Node */
    typedef boost::shared_ptr<Node> NodePtr;

    /** \brief Constructor */
    Node();
  
    /** \brief Destructor, virtual so that the children of derived objects gets destroyed properly */
    virtual ~Node();

    /** 
     * \brief Get a const pointer to the parent node
     * \return pointer to parent, or NULL if this node has no parent
     */
    Node* getParent() const;

    /**
     * \brief Get a pointer to a child
     * \param[in] idx Child index (should be less than getChildrenCount())
     * \return pointer to child, or NULL if idx is invalid
     */ 
    Node* getChildren(uint idx) const;

    /**
     * \brief Get the number of children that this node has
     * \return number of children
     */
    uint getChildrenCount() const;

    /**
     * \brief Add a child to this node. 
     * \note The parent node will take care of memory deallocation of children nodes
     * \param[in] child pointer to child node
     */
    void addChild(Node* child);
    
    /**
     * \brief Add a child to this node. 
     * \param[in] child node
     */
    void addChild(Node &child);

  protected:

    /** \brief children nodes */
    std::vector<NodePtr> children_;

    /** \brief parent node */
    Node* parent_;

  };



  /**
   * \class MurtyNode
   * \brief A node used in the Murty class for k-best linear assignment.
   * \author Keith Leung
   */
  class MurtyNode : public Node
  {
  public:

    /** \brief Smart pointer to a MurtyNode */
    typedef boost::shared_ptr<MurtyNode> Ptr;
    /** \brief A linear assignment */
    typedef boost::shared_array<int> Assignment;

    /** \brief Constructor 
     *  \param[in] id The partition id number for this node. Default is -1
     *  \param[in] n_children_exp the number of children this node is expected to have,
     *  used for more efficient memory allocation. Not mandatory.
     */
    MurtyNode(int id = -1, int n_children_exp = -1);

    /** \brief Destructor */
    ~MurtyNode();

    /** \brief Comparison operator for sorting */
    bool operator< (const MurtyNode& other);

    /** \brief Set the linear assignment represented by this node */
    void setAssignment(Assignment &a, double s);

    double getScore() const ;

    Assignment getAssignment() const;

    int getPartitionID() const;

  private:
  
    int id_; /**< \brief id number of this node */
    Assignment a_; /**< \brief Linear assignment */
    double s_; /**< \brief score of the linear assignment */

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

  /** \brief A linear assignment */
  typedef boost::shared_array<int> Assignment;
   
  /** \brief constructor 
   *  \param[in] C square score matrix from which we will find the best assignments.
   *  \param[in] n dimension of C
   *  \param[in] bigNum A number whose negative is used to replace entries in the cost matrix that cannot be used for assignment
   */ 
  Murty(double** C, unsigned int n, double bigNum = 10000);
  
  /** destructor */
  ~Murty();

  /**
   * \brief Get the best score.
   * \note The best score = 0 when object is first initiated.
   * \return the best score
   */
  double getBestScore() const;

  /** Find the next best assignment 
   *  \param[out] assignment resulting next best assignment
   *  \param[out] score the assignment score 
   *  \return k, the k-best solution index, or -1 if no more solutions are available
   */
  int findNextBest( Assignment &assignment, double &score);

  /** \brief Set the cost matrix block for real assignments. 
   *
   *  This is useful where we need to consider cases in whcih it is possible for elements
   *  to be unassigned. For example, in pairing landmarks with measurements, landmarks may be 
   *  misdetected, and measurements may be false alarm. To use this function, the origin (ideal)
   *  m x n assignment matrix should be augmented to be a (m+n) x (n+m) matrix, such that the
   *  ideal matrix is in the upper left. The upper right augmented block should represent miss-detections 
   *  and the lower left augmented block should represent false alarms. 
   *  \param[in] nR number of rows in the real assignment block
   *  \param[in] nC number of columns in the real assignment block
   */
  void setRealAssignmentBlock(unsigned int nR, unsigned int nC);

private:
  
  MurtyNode::Ptr root_; /**< The best linear assignment */
  int k_; /**< number of solutions found */
  double** C_; /**< the nxn score matrix */
  double** C_t_; /**< temporary score matrix */
  unsigned int n_; /**< dimension of C_ and C_t_*/
  HungarianMethod lam_; /**< Linear assignment method */
  std::priority_queue<MurtyNode*, std::vector<MurtyNode*>, MurtyNodeCompare> pq; /**< Priority queue for easy retrieval of node with highest score */
  double bestScore_;
  double bigNumber_;
  unsigned int realAssign_nC_; /**< dimension of ideal assignment block */
  unsigned int realAssign_nR_; /**< dimention of ideal assignment block */

};

}

#endif

