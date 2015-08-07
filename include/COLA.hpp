/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2015, Keith Leung, Felipe Inostroza
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

#ifndef COLA_HPP
#define COLA_HPP

#include "OSPA.hpp"

namespace rfs{

  /**
   * \class COLA
   * \f[ e_{\textrm{COLA}} = \frac{1}{c} \Bigg[ \min_{j \in \left\{ 1 \ldots \left|\mathcal{M}_g\right|\right\}} \sum_{i=1}^{\left|\mathcal{M}_k\right|} \min \left(d_{i,j},c\right)^p + c^p \Big| \left|\mathcal{M}_g\right| - \left|\mathcal{M}_k\right| \Big| \Bigg] ^{\frac{1}{p}} \f]
   * \brief A class for computing the Cardinalized Optimal Linear Assignment metric given 2 sets
   * \tparam SetElementType The class representing an element of the sets to be compared. The operator - must be defined in this class to return the distance between elements
   * \author Keith Leung
   */
  template<class SetElementType>
  class COLA : public OSPA<SetElementType>
  {

  public:
 
    /** Constructor 
     * \param set1 a std::vector containing the elements in set 1
     * \param set2 a std::vector containing the elements in set 2
     * \param cutoff value c
     * \param order value p
     */
    COLA(std::vector<SetElementType> &set1,
	 std::vector<SetElementType> &set2,
	 double cutoff,
	 double order);

    /** Destructor */
    ~COLA();

    /** Calculate the COLA error
     * \param[out] e_dist distance error component (before raising to power p_ ) 
     * \param[out] e_card cardinality error component (before raising to power p_ )
     * \param[in] report True to generate report to stdout
     * \return COLA error
     */
    double calcError(double *e_dist = NULL, double *e_card = NULL, bool report = false);

  };

  // IMPLEMENTATION
  
  template<class SetElementType>
  COLA<SetElementType>::COLA(std::vector<SetElementType> &set1,
			     std::vector<SetElementType> &set2,
			     double cutoff,
			     double order) : 
    OSPA<SetElementType>(set1, set2, cutoff, order)
  {
  }
    

  template<class SetElementType>
  COLA<SetElementType>::~COLA(){}

  template<class SetElementType>
  double COLA<SetElementType>::calcError(double *e_dist, double *e_card, bool report){

    double e_ospa = OSPA<SetElementType>::calcError(e_dist, e_card, report);
    double cost_ = e_ospa * pow(OSPA<SetElementType>::n_, 1.0/OSPA<SetElementType>::p_) / OSPA<SetElementType>::c_;

    return cost_;

  }


}

#endif
