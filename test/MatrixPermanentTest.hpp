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

#include <gtest/gtest.h>

#include "MatrixPermanent.hpp"

class MatrixPermanentTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  MatrixPermanentTest(){}
  
  /** Destructor */
  virtual ~MatrixPermanentTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

};

////////// Test Cases //////////

TEST_F(MatrixPermanentTest, permanentCalcTest){

  // Solutions according to: 
  // Nijenhuis, Albert, and Wilf, Herbert S.
  // Combinatorial Algorithms For Computers and Calculators 2nd Edition
  // Academic Press 1978
  
  // Matrices to be tested have zero diagonals and ones everywhere else

  double expected[13];
 
  expected[2] = 1;
  expected[3] = 2;
  expected[4] = 9;
  expected[5] = 44;
  expected[6] = 265;
  expected[7] = 1854;
  expected[8] = 14833;
  expected[9] = 133496;
  expected[10] = 1334961;
  expected[11] = 14684570;
  expected[12] = 176214841;

  for(int i = 2; i <= 12; i++){
    
    Eigen::MatrixXd A(i,i);
    A = Eigen::MatrixXd::Ones(i,i) - Eigen::MatrixXd::Identity(i,i);
    double actual = rfs::MatPerm::calc( A );
    ASSERT_DOUBLE_EQ(expected[i], actual);

  }

}
 
