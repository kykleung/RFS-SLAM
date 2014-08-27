#include <assert.h>
#include "CostMatrix.hpp"
#include <stdio.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace rfs
{

/// General Cost Matrix Implementation

CostMatrixGeneral::CostMatrixGeneral(double** &C, unsigned int nR, unsigned int nC) : nR_(nR),
										      nC_(nC),
										      components_i(NULL),
										      components_j(NULL),
										      components_mat(NULL),
										      components_rowIdx(NULL),
										      components_colIdx(NULL),
										      components_extendedAllocate_(NULL),
										      combinedZeroPartition_(-1),
										      nPartitions_(0)
{
  C_ = new double* [nR_];
  for(int i = 0; i < nR_; i++){
    C_[i] = new double [nC_];
  }
  C = C_;
}

CostMatrixGeneral::~CostMatrixGeneral(){

  for(int i = 0; i < nR_; i++){
    delete[] C_[i];
  }
  delete[] C_;

  if(components_mat != NULL){

    for(int p = 0; p < nPartitions_; p++){

      int nR = components_i[p].size();
      if(components_extendedAllocate_[p])
	nR += components_j[p].size();

      if((*components_mat)[p] != NULL){
	for(int r = 0; r < nR; r++){
	  delete[] (*components_mat)[p][r];
	}
	delete[] (*components_mat)[p];
      }
    }
    delete components_mat;
  }

  if(components_rowIdx != NULL){
    for(int p = 0; p < nPartitions_; p++){
      if((*components_rowIdx)[p] != NULL){
	delete[] (*components_rowIdx)[p];
      }
    }
    delete components_rowIdx;
  }

  if(components_colIdx != NULL){
    for(int p = 0; p < nPartitions_; p++){
      if((*components_colIdx)[p] != NULL){
	delete[] (*components_colIdx)[p];
      }
    }
    delete components_colIdx;
  }

  if(components_i != NULL){
    delete[] components_i;
  }
  if(components_j != NULL){
    delete[] components_j;
  }
  if(components_extendedAllocate_ != NULL){
    delete[] components_extendedAllocate_;
  }

}

void CostMatrixGeneral::getCostMatrix(double** &C, unsigned int &nR, unsigned int &nC){
  C = C_;
  nR = nR_;
  nC = nC_;
}


unsigned int CostMatrixGeneral::partition(){

  if(nPartitions_ != 0)
    return nPartitions_;

  // partitioning using connected-component analysis
  ::boost::adjacency_list< ::boost::vecS, ::boost::vecS, ::boost::undirectedS> G;
  for(int i = 0; i < nR_; i++){
    for(int j = 0; j < nC_; j++){
      if (C_[i][j] != 0 ){
	::boost::add_edge(i, j+nR_, G);
      }
    }
  }
  ::boost::add_edge(nR_+nC_-1, nR_+nC_-1, G); // add edge to self so the graph has the desired number of vetices

  std::vector<int> cc_results( nR_+nC_, -1 );
  int ncc = ::boost::connected_components(G, &cc_results[0]);

  if(components_i != NULL){
    delete[] components_i;
  }
  components_i = new std::vector<unsigned int> [ncc];
  for(int i = 0; i <  nR_; i++){
    components_i[cc_results[i]].push_back(i);
  }
  if(components_j != NULL){
    delete[] components_j;
  }
  components_j = new std::vector<unsigned int> [ncc];
  for(int j = nR_; j < nR_+nC_; j++){
    components_j[cc_results[j]].push_back(j - nR_);
  }

  combinedZeroPartition_ = -1;
  int nMergedZeroPartitions = 0; 
  for(int n = 0; n < ncc; n++){
    if(components_i[n].size() == 0 || components_j[n].size() == 0){
      
      if(combinedZeroPartition_ == -1){
	combinedZeroPartition_ = n;
      }else if(components_i[n].size() != 0){
	components_i[combinedZeroPartition_].push_back(components_i[n][0]);
	nMergedZeroPartitions++;
      }else{
	components_j[combinedZeroPartition_].push_back(components_j[n][0]);
	nMergedZeroPartitions++;
      }

    }
  }

  nPartitions_ = ncc - nMergedZeroPartitions;

  components_mat = new std::vector<double**>;
  components_rowIdx = new std::vector<unsigned int*>;
  components_colIdx = new std::vector<unsigned int*>;
  components_extendedAllocate_ = new bool [nPartitions_]; 
  
  components_mat->resize(nPartitions_, NULL);
  components_rowIdx->resize(nPartitions_, NULL);
  components_colIdx->resize(nPartitions_, NULL);

  return nPartitions_;

}


bool CostMatrixGeneral::getPartitionSize(int p, unsigned int &nRows, unsigned int &nCols){
  
  if( components_mat == NULL ){
    nRows = 0;
    nCols = 0;
  }
  nRows = components_i[p].size();
  nCols = components_j[p].size();
 
  if(p == combinedZeroPartition_)
    return false;
  return true;
 
}

bool CostMatrixGeneral::getPartition(int p, double** &Cp, unsigned int &nRows, unsigned int &nCols,
				     unsigned int* &row_indices, unsigned int* &col_indices, bool allocateExtendedMatrix){

  if( components_mat == NULL ){
    Cp = NULL;
    nRows = 0;
    nCols = 0;
    row_indices = NULL;
    col_indices = NULL;
    return false;
  }

  components_extendedAllocate_[p] = allocateExtendedMatrix;
  nRows = components_i[p].size();
  nCols = components_j[p].size();

  if((*components_mat)[p] != NULL){
    Cp = (*components_mat)[p];
    row_indices = (*components_rowIdx)[p];
    col_indices = (*components_colIdx)[p];
  }else{

    unsigned int nR = allocateExtendedMatrix ? nRows + nCols : nRows;
    unsigned int nC = allocateExtendedMatrix ? nRows + nCols : nCols;
    
    Cp = new double* [nR];
    for(int i = 0; i < nR; i++){
      Cp[i] = new double [nC];
      for(int j = 0; j < nCols; j++){
	if(i < nRows)
	  Cp[i][j] = C_[ components_i[p][i] ][ components_j[p][j]];
      }  
    }
    (*components_mat)[p] = Cp;
  
    row_indices = new unsigned int[nRows];
    for(int i = 0; i < nRows; i++){
      row_indices[i] = components_i[p][i];
    }
    (*components_rowIdx)[p] = row_indices;

    col_indices = new unsigned int[nCols];
    for(int j = 0; j < nCols; j++){
      col_indices[j] = components_j[p][j];
    }
    (*components_colIdx)[p] = col_indices;

  }
 
  if(p == combinedZeroPartition_)
    return false;
  return true;
}




/// Cost Matrix Implementation

CostMatrix::CostMatrix(double** &C, int nDim) : CostMatrixGeneral(C, nDim, nDim), 
					       n_(nDim), 
					       C_reduced_(NULL), 
					       n_reduced_(0),
					       reducedCostMatAvailable_(false),
					       i_remap_(NULL),
					       j_remap_(NULL)
					       
{

  a_fixed_.resize(nDim, -1);
  a_fixed_reverse_.resize(nDim, -1);
}

CostMatrix::~CostMatrix(){

  if( n_reduced_ > 0){
    for(int i = 0; i < n_reduced_; i++){
      delete[] C_reduced_[i];
    }
    delete[] C_reduced_;
  }
  if(i_remap_ != NULL)
    delete[] i_remap_;
  if(j_remap_ != NULL)
    delete[] j_remap_;

}

void CostMatrix::reduce(double lim, bool minVal){

  int nMatch_i[n_]; // number of possible matches
  int nMatch_j[n_];
  int best_score[n_];
  for(int n = 0; n < n_; n++){
    nMatch_i[n] = 0;
    nMatch_j[n] = 0;
    best_score[n] = lim;
  }

  for(int i = 0; i < n_; i++){
    for(int j = 0; j < n_; j++){

      if( (C_[i][j] >= lim && !minVal) || (C_[i][j] <= lim && minVal) ){ 
	C_[i][j] = lim;
      }else{
	nMatch_i[i]++;
	nMatch_j[j]++;

	if( nMatch_i[i] == 1 && nMatch_j[j] == 1){
	  a_fixed_[i] = j;
	  a_fixed_reverse_[j] = i;
	}

	if( nMatch_i[i] > 1){
	  a_fixed_[i] = -1;
	}
	if( nMatch_j[j] > 1){
	  a_fixed_reverse_[j] = -1;
	}
      }
      
    }
  }

  for(int n = 0; n < n_; n++){
    if( nMatch_j[ a_fixed_[n] ] != 1 ){
      a_fixed_[n] = -1;
    }
    if( nMatch_i[ a_fixed_reverse_[n] ] != 1){
      a_fixed_reverse_[n] = -1;
    }
    if( a_fixed_[n] == -1){
      i_reduced_.push_back(n);
    }
    if( a_fixed_reverse_[n] == -1){
      j_reduced_.push_back(n);
    }
  }

  assert( i_reduced_.size() == j_reduced_.size());
  n_reduced_ = i_reduced_.size();

  if(n_reduced_ != 0){

    C_reduced_ = new double* [n_reduced_];
    for(int i = 0; i < n_reduced_; i++ ){
      C_reduced_[i] = new double [n_reduced_];
      for(int j = 0; j < n_reduced_; j++ ){
	C_reduced_[i][j] = C_[ i_reduced_[i] ][ j_reduced_[j] ];
      } 
    }

    if(n_reduced_ == 1){

      a_fixed_[ i_reduced_[0] ] = j_reduced_[0];
      n_reduced_ = 0;

    }
  }

  reducedCostMatAvailable_ = true;
}

int CostMatrix::getCostMatrixReduced(double** &C, int* &fixedAssignments, double* fixedScore, int* &i_remap, int* &j_remap ){

  if(!reducedCostMatAvailable_)
    return -1;

  if(i_remap_ == NULL)
    i_remap_ = new int[n_];
  if(j_remap_ == NULL) 
    j_remap_ = new int[n_];

  C = C_reduced_;
  *fixedScore = 0;
  for(int n = 0; n < n_; n++){
    fixedAssignments[n] = a_fixed_[n];
    if(fixedAssignments[n] != -1){
      *fixedScore += C_[n][fixedAssignments[n]];
    }
  }

  for(int n = 0; n < n_reduced_; n++){
    i_remap_[n] = i_reduced_[n];
    j_remap_[n] = j_reduced_[n];
  }
  i_remap = i_remap_;
  j_remap = j_remap_;

  return n_reduced_;
}

}
