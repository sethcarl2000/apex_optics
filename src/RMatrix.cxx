//*-- Author :    Seth Hall;  22 Jul 2024

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMatrix
//
// This is a very bare-bones Lin. algebra class, meant to deal primarily
// with square-matricies. Specifically, it's designed to work with
// the ROOT::vecd framework. 
//
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include "RMatrix.h"

using namespace std;

using vecd = std::vector<double>; 

//_______________________________________________________________________________
//_______________________________________________________________________________
RMatrix::RMatrix(unsigned int nr, unsigned int nc, double init)
  : fnCols(nc),
    fnRows(nr),
    f_isSquare(nr==nc), 
    f_reportSingular(true),
    f_n_elems(nc*nr)
{
  fElems = vector<double>(f_n_elems, init); 
}
//_______________________________________________________________________________
RMatrix::RMatrix(unsigned int nr, unsigned int nc, const vecd &array)
  : fnCols(nc),
    fnRows(nr),
    f_isSquare(nr==nc), 
    f_reportSingular(true),
    f_n_elems(nc*nr)
{
  if (f_n_elems != array.size()) {
    fprintf(stderr, "Warning in <RMatrix::RMatrix(const vecd &)>: "  
	  "Array size does not match given matrix dims! Initialzed as all zeros.\n");
    fElems = vector<double>(f_n_elems, 0.); 
  } else { 
    fElems = array;
  }
}
//_______________________________________________________________________________
RMatrix::RMatrix(const unsigned int nr,const unsigned int nc, const vecd *array)
  : fnCols(nc),
    fnRows(nr),
    f_isSquare(nr==nc), 
    f_reportSingular(true),
    f_n_elems(nc*nr)
{
  if (nr*nc != array->size()) {
    fprintf(stderr, "Warning in <RMatrix::RMatrix(const vecd*)>: "  
	  "Array size does not match given matrix dims! Initialzed as all zeros.\n");
    fElems = vector<double>(f_n_elems, 0.); 
  } else { 
    fElems = *array;
  }
}
//_______________________________________________________________________________
RMatrix::~RMatrix() 
{
  fElems.clear(); 
}
//_______________________________________________________________________________
/*RMatrix::RMatrix(const RMatrix &mat)
  : TObject(),
    fnCols(mat.GetNCols()),
    fnRows(mat.GetNRows()),
    f_isSquare(fnCols==fnRows),
    f_reportSingular(true),
    fElems(mat.Data())
{
  //copy-constructor
  } */
//_______________________________________________________________________________
double& RMatrix::at(unsigned int i, unsigned int j)
{
  if (i>=GetNRows() || j>=GetNCols()) {
    fprintf(stderr, "Error in <RMatrix::at>: "
    "Invalid element access attempt: (%i,%i), last elem is (%i,%i)\n",
	  i,j, GetNRows()-1,GetNCols()-1);
    return *(new double(std::nan("1")));
  }
  return fElems[GetNCols()*i + j]; 
}
//_______________________________________________________________________________
double RMatrix::at(unsigned int i, unsigned int j) const
{
  if (i>=GetNRows() || j>=GetNCols()) {
    fprintf(stderr, "Error in <RMatrix::at>: "
    "Invalid element access attempt: (%i,%i), last elem is (%i,%i)\n",
	  i,j, GetNRows()-1,GetNCols()-1);
    return *(new double(std::nan("1")));
  }
  return fElems[GetNCols()*i + j]; 
}
//_______________________________________________________________________________
vecd RMatrix::operator*(const vecd &rhs) const
{
  //check vector size
  if (rhs.size() != GetNCols()) {
    fprintf(stderr, "Error in <RMatrix::operator*>: "
    "Tried to multiply a column vector of size %i by a matrix of size %i X %i.\n",
    (int)rhs.size(), GetNRows(), GetNCols() );
    return {};
  }
  
  vecd out(fnRows,0.);
  
  for (unsigned int i=0; i<GetNRows(); i++) {
    for (unsigned int j=0; j<GetNCols(); j++) { out[i] += fElems[GetNCols()*i + j] * rhs[j]; }
  }
  
  return out; 
} 
//_______________________________________________________________________________
RMatrix RMatrix::operator+(const RMatrix &rhs) const
{
  //check matrix sizes
  if (rhs.GetNCols() != GetNCols() ||
      rhs.GetNRows() != GetNRows()) {
    fprintf(stderr, "Error in <RMatrix::operator+>: "
    "Tried to add matrix of size %i X %i by a matrix of size %i X %i.\n",
    (int)GetNRows(),    (int)GetNCols(),
	  (int)rhs.GetNRows(),(int)rhs.GetNCols());
    return RMatrix();
  }

  //make a copy of this matrix
  RMatrix sum = *this;
  
  for (unsigned int i=0; i<GetNRows(); i++) 
    for (unsigned int j=0; j<GetNCols(); j++) sum.get(i,j) += rhs.at(i,j);
  
  return sum; 
} 
//_______________________________________________________________________________
vecd RMatrix::Solve(const vecd &B) const
{
  //check if we're trying to call this when this matrix is not square
  if (!f_isSquare) {
    fprintf(stderr, "Error in <RMatrix::Solve>: "
    "Tried to invert a non-square matrix!"); 
    return {};
  }

  const unsigned int N = GetNRows(); 
  
  //takes the NxN matrix A, and the N-vector B, and solves it. 
  //performs numerical LU factorization 
  
  //ASSUMES MATRIX IS NONSINGULAR!!!!
  
  //check vector size
  if (B.size() != N) {
    fprintf(stderr, "Error in <RMatrix::Solve>: "
    "Tried to Solve a lin. system with a column vector of size %i and a matrix of size %i X %i.\n",
    (int)B.size(), GetNRows(), GetNCols() );
    return {};
  }

  //the U-matrix starts as a copy of the 'A' input-matrix
  RMatrix U = *this;  
  //we initialize the L-matrix with all zeros
  RMatrix L(N,N, 0.);  
  
  //Create the U (upper triangular) and L (lower triangular) matrices
  for (unsigned int ii=0; ii<N; ii++) { 
    //loop over a (successivley smaller) sub-matrix
    double a_00 = U.at(ii,ii); 
    for (unsigned int i=ii+1; i<N; i++) { 
      double a_i0 = U.at(i,ii); 
      L.at(i,ii) = a_i0/a_00; 
      for (unsigned int j=ii; j<N; j++) { 
        U.at(i,j) += ( -a_i0/a_00 ) * U.at(ii,j); 
      }
    }
  }
  double det = 1.; 
  for (unsigned int ii=0; ii<N; ii++) { L.at(ii,ii)=1.; det *= U.at(ii,ii); }

  //det is NaN 
  if ( det != det ) {
    if (f_reportSingular) fprintf(stderr, "Error in <RMatrix::Solve>: Determinant is NAN."); 
    return {}; 
  }
  
  vecd y(N,0.);

  //solve the system Ly = b
  for (unsigned int i=0; i<N; i++) { 
    y[i] = B[i]; 
    for (unsigned int j=0; j<i; j++) y[i] += -L.at(i,j)*y[j];
  }

  vecd X(N,0.); 
  
  //now solve Ux = y
  for (int i=N-1; i>=0; i--) { 
    X[i] = y[i]; 
    for (int j=N-1; j>i; j--) { 
      X[i]  +=  - U.at(i,j) * X[j]; 
    }
    X[i] *= 1./U.at(i,i); 
  }
  return X; 
}
//_______________________________________________________________________________
//_______________________________________________________________________________
//_______________________________________________________________________________
void RMatrix::Print()
{
  //print the elems of this matrix (for debugging purposes)
  printf("Matrix: %ix%i\n",GetNRows(),GetNCols()); 

  for (unsigned int i=0; i<GetNRows(); i++) {
    for (unsigned int j=0; j<GetNCols(); j++) 
      printf(" %+0.3e", fElems[i*GetNCols() + j] );
    cout << endl;
  } 
}
//_______________________________________________________________________________

ClassImp(RMatrix)
