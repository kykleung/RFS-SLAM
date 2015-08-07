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

#include "SpatialIndexTree.hpp"
#include "Landmark.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>


class SpatialIndexTreeTest : public ::testing::Test{

protected:

  typedef rfs::Landmark<3> Pt;
  typedef boost::shared_ptr<Pt> Ptr;

  /** Constructor for setting up each test */
  SpatialIndexTreeTest(){

    coord_max = 1000;
    coord_min = -coord_max;
    box_size_min = 2;
    tree = rfs::SpatialIndexTree<3, rfs::Landmark<3> >(box_size_min);
    nLandmarks = 5000;
  }
  
  /** Destructor */
  virtual ~SpatialIndexTreeTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //
  rfs::SpatialIndexTree<3, rfs::Landmark<3> > tree;
  double coord_max;
  double coord_min;
  double box_size_min;
  std::vector<Ptr> map_added;
  std::vector<Ptr> map_removed;
  int nLandmarks;

};

////////// Test Cases //////////

TEST_F(SpatialIndexTreeTest, simplePopultateTreeTest){
 
  Ptr m1(new Pt);
  Ptr m2(new Pt);
  Ptr m3(new Pt);

  (*m1)[0] = 330.548;
  (*m1)[1] = 321.398;

  (*m2)[0] = -35.452;
  (*m2)[1] = 322.292;

  (*m3)[0] = 746.554;
  (*m3)[1] = 71.070;

  map_added.push_back(m1);
  tree.addData(m1);
  //tree.exportASCII("tree1.dat");
  ASSERT_EQ(1, map_added.size());
  ASSERT_EQ(1, tree.getDataSize() );

  map_added.push_back(m2);
  tree.addData(m2);
  //tree.exportASCII("tree2.dat");
  ASSERT_EQ(2, map_added.size());
  ASSERT_EQ(2, tree.getDataSize() );

  map_added[1] = map_added.back();
  map_added.pop_back();
  tree.removeData(m2);
  //tree.exportASCII("tree3.dat");
  ASSERT_EQ(1, map_added.size());
  ASSERT_EQ(1, tree.getDataSize() );

  map_added.push_back(m3);
  tree.addData(m3);
  //tree.exportASCII("tree4.dat");
  ASSERT_EQ(2, map_added.size());
  ASSERT_EQ(2, tree.getDataSize() );

}

TEST_F(SpatialIndexTreeTest, populateTreeTest){

  boost::variate_generator< ::boost::mt19937, ::boost::uniform_real<double> > 
    gen(boost::mt19937(time(NULL)), ::boost::uniform_real<double>(coord_min, coord_max));

  // Randomly populate map space
  for(int k = 0; k < 100; k++){
    //std::cout << "================================\n";
    while(map_added.size() < nLandmarks){
      
      bool add = true;
      if(map_added.size() > 0 && gen() > (coord_max - coord_min) * 0.55 + coord_min ){
	add = false;
      }

      if(add){
	Ptr m = Ptr(new Pt);
	for(int i = 0; i < m->getNDim(); i++){
	  (*m)[i] = gen();
	}
	map_added.push_back(m);
	// std::cout << "add " << m->get().transpose() << "\n";
	tree.addData(m);
      }else{
	boost::variate_generator< ::boost::mt19937, boost::uniform_int<> > 
	  gen_int(::boost::mt19937(rand()), ::boost::uniform_int<>(0, map_added.size() - 1));
	int remove_idx = gen_int();
	// std::cout << "remove " <<  remove_idx << "/" << map_added.size() << " [" << map_added[remove_idx]->get().transpose() <<  "]\n"; 
	map_removed.push_back( map_added[remove_idx] );
	map_added[remove_idx] = map_added.back();
	map_added.pop_back();
	tree.removeData(map_removed.back());
      }

      ASSERT_EQ(map_added.size(), tree.getDataSize());
    }
    ASSERT_EQ(map_added.size(), tree.getDataSize());

    map_added.clear();
    map_removed.clear();

    //tree.exportASCII("tree.dat");
    tree.clear();
    ASSERT_EQ(0, tree.getDataSize());
  }
 
}


TEST_F(SpatialIndexTreeTest, closestNeighbourSearchTest){

  boost::variate_generator< ::boost::mt19937, ::boost::uniform_real<double> > 
    gen(boost::mt19937(time(NULL)), ::boost::uniform_real<double>(coord_min, coord_max));

  // Randomly populate map space
  while(map_added.size() < nLandmarks){
  
    Ptr m = Ptr(new Pt);
    for(int i = 0; i < m->getNDim(); i++){
      (*m)[i] = gen();
    }
    map_added.push_back(m);
    //std::cout << "add " << m->get().transpose() << "\n";
    tree.addData(m);

  }
  
  Ptr q = Ptr(new Pt); // query point
  for(int k = 0; k < 1000; k++ ){
    for(int i = 0; i < q->getNDim(); i++){
      (*q)[i] = gen();
    }
    
    std::vector<Ptr> m_closest_bruteForce;
    m_closest_bruteForce.push_back(map_added[0]);
    double dist = (map_added[0]->get() - q->get()).norm();
    for(int i = 1; i < map_added.size(); i++){
      double d_temp = (map_added[i]->get() - q->get()).norm();
      if(d_temp < dist){
	m_closest_bruteForce.clear();
	m_closest_bruteForce.push_back(map_added[i]);
	dist = d_temp;
      }else if(d_temp == dist){
	m_closest_bruteForce.push_back(map_added[i]);
      }
    }

    Ptr m_closest_tree;
    m_closest_tree = tree.queryClosestPt(q->get());

    bool closestFound = false;
    for(int i = 0; i < m_closest_bruteForce.size(); i++){
      if(m_closest_bruteForce[i] == m_closest_tree)
	closestFound = true;
    }
    ASSERT_TRUE(closestFound);
  }


}

TEST_F(SpatialIndexTreeTest, boxSearchTest){

  boost::variate_generator< ::boost::mt19937, ::boost::uniform_real<double> > 
    gen(boost::mt19937(time(NULL)), ::boost::uniform_real<double>(coord_min, coord_max));
  
  // Randomly populate map space
  while(map_added.size() < nLandmarks){
  
    Ptr m = Ptr(new Pt);
    for(int i = 0; i < m->getNDim(); i++){
      (*m)[i] = gen();
    }
    map_added.push_back(m);
    //std::cout << "add " << m->get().transpose() << "\n";
    tree.addData(m); 
  }
  ASSERT_EQ(map_added.size(), tree.getDataSize());
  
  Ptr q_max = Ptr(new Pt); // query box max corner
  Ptr q_min = Ptr(new Pt); // query box min corner
  for(int k = 0; k < 1000; k++ ){
    
    // randomly generate query box
    for(int i = 0; i < q_max->getNDim() ; i++){
      do{
	  (*q_max)[i] = gen();
	  (*q_min)[i] = gen();
      }while((*q_max)[i] <= (*q_min)[i]);
    }    
    // std::cout << "query box: " << q_min->get().transpose() << " " << q_max->get().transpose() << std::endl;
	  
    std::vector<Ptr> query_result_bruteForce;
    for(int i = 0; i < map_added.size(); i++){
      if( (q_max->get() - map_added[i]->get()).minCoeff() >= 0 && (q_min->get() - map_added[i]->get()).maxCoeff() <= 0 ){
	query_result_bruteForce.push_back(map_added[i]);
      }
    }
    
    std::vector<Ptr> query_result_tree;
    tree.queryDataInBox(q_min->get(), q_max->get(), query_result_tree);

    // sizes of the results by the two methods should be identical
    EXPECT_EQ(query_result_bruteForce.size(), query_result_tree.size());

    // content of the results should be identical, although ordering may be different
    for(int i = 0; i < query_result_bruteForce.size(); i++){
      bool found = false;
      for(int j = 0; j < query_result_tree.size(); j++){
	if( query_result_bruteForce[i] == query_result_tree[j] ){
	  query_result_tree[j] = query_result_tree.back();
	  query_result_tree.pop_back();
	  found = true;
	  break;
	}
      }
      //if(!found) std::cout << "No match for " << query_result_bruteForce[i]->get().transpose() << std::endl;
      ASSERT_TRUE(found);
    }
    ASSERT_EQ(0, query_result_tree.size());
    
  }
}

