#ifndef RMatrix_h_
#define RMatrix_h_

//////////////////////////////////////////////////////////////////////////
//
// RMatrix
//
// A jankly little lin-algebra class which is able to work (implicitly)
// with the ROOT::RVec<double> class, 'cause they're a lot more convenient
// to work with than TVectorD's are. 
//
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <memory> 
#include <vector>
#include <ROOT/RVec.hxx>

#include "TROOT.h"
#include "TObject.h"


//________________________________________________________________________________
class RMatrix : public TObject {
public:
  
  RMatrix(unsigned int nr=1, unsigned int nc=1, double init=0.); 

  RMatrix(unsigned int nr, unsigned int nc, const ROOT::RVec<double> &array);
  
  RMatrix(unsigned int nr, unsigned int nc, const ROOT::RVec<double> *array);

  //copy constructor
  RMatrix(const RMatrix& rhs) noexcept; 
  
  //move constructor
  RMatrix(RMatrix&& rhs) noexcept; 
  
  //move assignment operator
  RMatrix& operator=(RMatrix&& rhs) noexcept; 
  
  //copy constructor
  //RMatrixD(const RMatrixD &mat); 
  
  static RMatrix OuterProduct(const ROOT::RVec<double> &u, const ROOT::RVec<double> &v); 

  ~RMatrix();  

  ROOT::RVec<double> Solve(const ROOT::RVec<double> &B) const;

  double Determinant() const; 

  //operators: 
  //multiplication by ROOT::RVec<double>
  ROOT::RVec<double> operator*(const ROOT::RVec<double> &rhs) const; 

  //multiplication by scalar 
  RMatrix            operator*(double mult) const; 

  //addition operator (addition of two RMatirx objects, to be stored in a third)
  RMatrix            operator+(const RMatrix &rhs) const; 

  //add a matrix to this one. this is optimized for speed. 
  void               operator+=(RMatrix &rhs);       

  //perform matrix multiplication
  RMatrix operator*(RMatrix& rhs); 



  
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
  inline bool ReportSingular() const { return f_reportSingular; }
  inline void Set_report_singular(bool _val) { f_reportSingular=_val; }

  //return a copy to the data
  ROOT::RVec<double>& Data() { return fElems; };
   
  //return a const ptr to the underlying data 
  const ROOT::RVec<double> Data_cpy() const { return fElems; }

  //return a sqaure identity matrix
  static RMatrix Square_identity(int size); 

  
private:

  unsigned int fnCols, fnRows; 
  bool f_isSquare; 
  
  // (default==true) 
  // controls wheter or not an error message is printed when a singular matrix is encountered
  bool f_reportSingular; 
  
  int f_n_elems; 
  ROOT::RVec<double> fElems; 

  ClassDef(RMatrix,1)
}; 
//________________________________________________________________________________

#endif
