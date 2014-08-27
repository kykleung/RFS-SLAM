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

#ifndef LINEAR_ASSIGNMENT_HPP
#define LINEAR_ASSIGNMENT_HPP

#include "CostMatrix.hpp"
#include <queue>

namespace rfs{

/**
 * \class BruteForceLinearAssignment
 * This function finds all the linear assignments and orders them from best to worst.
 * It is intended for testing and checking the Hungarian method and Murty's k-best algorithm.
 * \brief Brute-force linear assignment 
 */
class BruteForceLinearAssignment
{

public:

  /** Constructor */
  BruteForceLinearAssignment();

  /** Destructor */
  ~BruteForceLinearAssignment();
  
  /**
   *  Run brute-force linear assignment by finding the cost of all possible linear assignments
   *  in a lexicographical order, and then returning the results in an ordered format with 
   *  respect to the cost / score.
   *  \param[in] C square cost / score matrix
   *  \param[in] n dimension of C
   *  \param[out] a ordered assignments a[k][n], where k is the assignment order
   *  \param[out] s ordered assignment scores
   *  \param[in] maxToMin ordering of results
   *  \return number of assignments
   */
  unsigned int run(double** C, int n, unsigned int** &a, double* &s, bool maxToMin = true);

private:

  /** \brief A linear assignment */
  struct assignment{
    unsigned int* a;
    double score;
    
    bool operator<(const assignment& rhs) const{
      if(score < rhs.score)
	return true;
      return false;
    }

  };

  std::priority_queue<assignment> pq_;
  unsigned int nAssignments_;
  unsigned int** a_;
  double* s_;
};

}

#endif
