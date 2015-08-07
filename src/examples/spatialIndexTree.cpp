#include <iostream>
#include <vector>

#include "SpatialIndexTree.hpp"

int main(int argc, char** argv){

  rfs::SpatialIndexTree<2> t(2);
  boost::shared_ptr< rfs::RandomVec<2> > px;
  rfs::RandomVec<2>::Vec xv;
  
  xv << 1, 1;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);

  xv << 3, 3;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);
  
  xv << 1.5, 1.5;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);
  
  xv << -1.5, -1.5;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);

  xv << -0.5, -0.5;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);
  
  xv << 5, -1;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);
  
  xv << -3, -1;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);
  
  std::cout << "remove data: " << t.removeData(px) << std::endl;

  xv << 4, -6;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);

  xv << -3, -3;
  px = boost::shared_ptr< rfs::RandomVec<2> >(new rfs::RandomVec<2>);
  px->set(xv);
  t.addData(px);

  t.exportASCII("tree.dat");

  // find closest point to
  xv << -10, 0;
  boost::shared_ptr< rfs::RandomVec<2> > qx = t.queryClosestPt(xv);
  std::cout << "closest point to\n" << xv << "\nis\n" << qx->get() << "\n";

  // find points in box
  std::vector< boost::shared_ptr< rfs::RandomVec<2>  >  > query_result;
  rfs::RandomVec<2>::Vec box_min;
  rfs::RandomVec<2>::Vec box_max;
  box_min << 2,-6;
  box_max << 6,-0.5;
  t.queryDataInBox(box_min, box_max, query_result);
  std::cout << "Query box \n" << box_min << "\n" << box_max << "\n\n";
  std::cout << "Query results:\n";
  for(int i = 0; i < query_result.size(); i++){
    std::cout << query_result[i]->get() << std::endl << std::endl;
  }
  
  return 0;
} 
