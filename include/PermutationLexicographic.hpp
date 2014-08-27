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

#ifndef PERMUTATION_LEXICOGRAPHIC
#define PERMUTATION_LEXICOGRAPHIC

namespace rfs{

/**
 * \class PermutationLexicographic
 * A class for generating assingments based on lexicographic ordering with or without clutter.
 * This class is mainly intended for assignments of measurements to landmarks, where measurements
 * can be outliers (i.e., non-associated with a landmark), and landmarks can be miss-detected
 * (i.e., not associated with a measurement)
 * \brief Generate lexicographical ordering
 */
class PermutationLexicographic
{
public:
  
  /**
   * Constructor 
   * \param[in] nM number of landmarks
   * \param[in] nZ number of measurements
   * \param[in] includeClutter if true, considers linear assignments with clutter and miss-detections. 
   * This will default to true if nM and nZ are not the same
   */
  PermutationLexicographic(unsigned int nM, unsigned int nZ, bool includeClutter = true);

  /**
   * Destructor 
   */
  ~PermutationLexicographic();

  /**
   * Find the next permutation in the lexicographic ordering
   * \param[out] permutation pointer to an array for reading the current permutation. 
   * Memory should be allocated by the caller. The size of the array needs to be nM if includeClutter = false,
   * and nM + nZ if includeClutter is true.
   * \return permutation number or 0 if there are no permutations remaining
   */
  unsigned int next(unsigned int* permutation);

private:

  unsigned int* o_; /**< ordering */
  unsigned int oSize_; /**< size of array o_ */
  unsigned int nM_; /**< number of landmarks */
  unsigned int nZ_; /**< number of measurements */
  unsigned int nP_; /**< number of permutations */
  bool last_; /**< last permutation reached */
};

}

#endif
