#ifndef NPoly_h_
#define NPoly_h_

//////////////////////////////////////////////////////////////////////////
//
//  NPoly
// 
//  A class which is a N-dimensional polynomal, with N-inputs and which
//  generates 1 output. Also has the ability to auto-generate elements
//  given a number of DoF and max polynomial 'order', via the 
//  AutoConstructPoly method. 
//
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TObject.h"
#include <ROOT/RVec.hxx> 
#include "RMatrix.h"


class NPoly : public TObject {
  
 public:

  NPoly() {}; 
  NPoly(int nDoF, int order=0); 
  ~NPoly();
  
  struct NPolyElem {
    ROOT::RVec<int>  powers; //powers for each (x,y,th,ph) in target-coordinates  
    double coeff;            //in most cases, coefficients are not needed per poly.

    //only the 'exponent signature' of these elements is compared; coefficient is not considered. 
    bool operator==(const NPolyElem& rhs) const { 
        const int my_size(powers.size()); 
        if (my_size != rhs.powers.size()) return false; 
        for (int i=0; i<my_size; i++) if (powers[i] != rhs.powers[i]) return false; 
        return true; 
    }
  };

  //return the value of each element evaluated against the input vector
  ROOT::RVec<double> Eval_noCoeff(const ROOT::RVec<double> &X) const;
  
  double Eval(const ROOT::RVec<double> &coeff, const ROOT::RVec<double> &X) const; 
  double Eval(const ROOT::RVec<double> &X) const;
  
  //compute the gradient w/r/t each of the input coordinates
  ROOT::RVec<double> Gradient(const ROOT::RVec<double> &coeff, const ROOT::RVec<double> &X) const; 
  
  //compute the gradient w/r/t each of the input coordinates
  ROOT::RVec<double> Gradient(const ROOT::RVec<double> &X) const; 

  //compute the hessian matrix (H_ij = [d^2 / dx_i dx_j]Poly) w/r/t the input coordinates
  RMatrix Hessian(const ROOT::RVec<double> &X) const;

  //handle elements + get information
  //Add element to the poly.
  void Add_element(const ROOT::RVec<int> &pows, double coefficient=1.);
  //Add element to the poly. 
  void Add_element(const NPoly::NPolyElem& elem); 
  
  //get the number of elements currently in the polynomial
  inline unsigned int Get_nElems() const { return fElems.size(); }
  
  //get this polynomial's input DoF
  inline unsigned int Get_nDoF() const { return fnDoF; }

  //Maximum exponent power from any element
  inline unsigned int Get_maxPower() const { return f_maxPower; }

 
  //Get a pointer to an individual element
  const NPolyElem*    Get_elem      (unsigned int i) const; 
  NPolyElem*          Get_elem      (unsigned int i); 
  //Get the RVec<int> of the powers of each input var for this element
  ROOT::RVec<int>     Get_elemPowers(unsigned int i) const;
  //get the coefficient of this element
  double              Get_elemCoeff (unsigned int i) const;

  NPolyElem*          Find_element(const ROOT::RVec<int>& powers); 

  //for debug purposes, print
  void Print() const;
  

  //these following methods will be used for symbolic computation; 
  // for example, symbolically computing the result of feeding the output of one NPolyArray into the input of another. 
  // I am experimenting with constructing optics models which are composed of several NPolyArray's 'chained together' 
  // in which the output of each NPolyArray in the chain is fed into the input of the next. 

  //symbolically compute the result of multiplying two NPoly's together.  
  NPoly operator*(const NPoly& rhs) const;  

  //symbolically compute the result of multiplying two NPoly's together.  
  NPoly operator+(const NPoly& rhs) const; 

  //this multiplies the coefficient of each element by a const number. 
  void operator*=(double val); 

  //raise a polynomial to an integer power
  static NPoly Pow(const NPoly& pol, const int pow); 
  
private:
  
  int fOrder,fnDoF;

  int f_maxPower; 
  
  ROOT::RVec<NPolyElem> fElems;
  
  //construct all possible polynomial elemts for which the sum of the exponents for all DoF is <= max_power. 
  void AutoConstructPoly(const int max_power, const int nDoF); 
  
  ClassDef(NPoly,1);   
};

#endif
