#ifndef RMatrix_h_
#define RMatrix_h_

//////////////////////////////////////////////////////////////////////////
//
// RMatrix
//
// A jankly little lin-algebra class which is able to work (implicitly)
// with the std::vector<double> class, 'cause they're a lot more convenient
// to work with than TVectorD's are. 
//
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <memory> 
#include <vector>

#include "TROOT.h"
#include "TObject.h"


//________________________________________________________________________________
class RMatrix : public TObject {
public:
  
  RMatrix(unsigned int nr=1, unsigned int nc=1, double init=0.); 

  RMatrix(unsigned int nr, unsigned int nc, const std::vector<double> &array);
  
  RMatrix(unsigned int nr, unsigned int nc, const std::vector<double> *array);
  
  //copy constructor
  //RMatrixD(const RMatrixD &mat); 
  
  ~RMatrix();  

  std::vector<double> Solve(const std::vector<double> &B) const;

  
  //multiplication by std::vector<double>
  std::vector<double> operator*(const std::vector<double> &rhs) const; 

  //adding two matrices
  RMatrix             operator+(const RMatrix &rhs) const; 
  
  inline unsigned int GetNCols() const { return fnCols; }
  inline unsigned int GetNRows() const { return fnRows; }

  //element-wise access
  double& at(unsigned int i,unsigned int j);
  double  at(unsigned int i,unsigned int j) const; 
    
  //same as at(i,j), but without bounds-checking (dangerous!)
  inline double &get(unsigned int i, unsigned int j) 
  { return fElems[GetNCols()*i + j]; }

  //print out elements
  void Print(); 
  
  //toggle whether or not this matrix will spit out an error if its singular 
  inline bool &ReportSingular() {return f_reportSingular;}

  //return a copy to the data
  std::vector<double> *Data() { return &fElems; }; 
  
private:

  unsigned int fnCols, fnRows; 
  bool f_isSquare; 
  
  // (default==true) 
  // controls wheter or not an error message is printed when a singular matrix is encountered
  bool f_reportSingular; 
  
  int f_n_elems; 
  std::vector<double> fElems; 

  ClassDef(RMatrix,1)
}; 
//________________________________________________________________________________

#endif
