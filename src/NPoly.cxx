//*-- Author :    Seth Hall   09-Jul-24

//////////////////////////////////////////////////////////////////////////
//     
// NPoly
// This is the backbone of the 'reverse optics' model.
// effectivley, all it is is a polynomial, with a appendable list of elements,
// defined on some R^n space.
// In the case of the reverse-optics, it will be defined on an R^4 space,
//  + using 3 polynomials in conjunction lets you make a f:R^4->R^3 map.
//  
//     
//////////////////////////////////////////////////////////////////////////

#include "NPoly.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iostream> 

using namespace ROOT::VecOps; 
using vecd = RVec<double>; 

using namespace std;

//_____________________________________________________________________________
NPoly::NPoly(const int nDoF, const int order) :
  fnDoF(nDoF), f_maxPower(0)
{  
  if (order>0) {
    fOrder = order;
    AutoConstructPoly(order,nDoF);
  }
}
//_____________________________________________________________________________
NPoly::~NPoly() {/*destructor*/};
//_____________________________________________________________________________
vecd NPoly::Eval_noCoeff(const vecd &X) const
{
  //check to make sure that the input vector is the right size
  if ((int)X.size() != Get_nDoF()) {
    Error("NPoly::Eval_noCoeff",
	  "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), Get_nDoF());
    return {};
  }

  vecd vals; vals.reserve(Get_nElems());
  
  const int max_pow = Get_maxPower(); 
  double X_pows[Get_nDoF()][max_pow + 1]; 
  for (int d=0; d<Get_nDoF(); d++) {
    for (int p=0; p<=max_pow; p++) X_pows[d][p] = pow(X[d], p); 
  }

  for (const NPolyElem &elem : fElems) {

    double elem_val=1.;

    for (int d=0; d<Get_nDoF(); d++) {
      if (elem.powers[d] > 0) elem_val *= X_pows[d][elem.powers[d]];
    }
    //check to see if this is the one element for which all pows. are zero
    vals.push_back( elem_val );
  }

  return vals;
}
//_____________________________________________________________________________
double NPoly::Eval(const vecd &coeff, const vecd &X) const
{
  if ((int)coeff.size() != Get_nElems()) { 
    Error("Eval(vecd,vecd)",
	  "Wrong number of coeffs. given: expected %i, recieved %i.",
	  Get_nElems(),
	  (int)coeff.size() );
    return -1e30;
  }
        
  if ((int)X.size() != Get_nDoF()) { 
    Error("Eval(vecd,vecd)", "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), Get_nDoF());
    return -1e30;
  }
      
  //now, actually evaluate
  double val=0.; 

  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff.
  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff.
  int i=0; 
  for (const auto &elem : fElems) {
    
    double elem_val=1.;
    
    for (int d=0; d<Get_nDoF(); d++) {
      if (elem.powers[d]>0) elem_val *= pow( X[d], elem.powers[d] );
    }
    //check to see if this is the one element for which all pows. are zero
    val += elem_val * coeff[i++]; 
	
  }//for (const auto &elem : fElems)

  /*for (unsigned int i=0; i<Get_nElems(); i++) { const auto& elem = fElems[i]; 
    
    double elem_val=1.;
    
    for (int d=0; d<Get_nDoF(); d++) {
      if (elem.powers[d]>0) elem_val *= pow( X[d], elem.powers[d] );
    }
    //check to see if this is the one element for which all pows. are zero
    val += elem_val * coeff[i]; 
	
  }//for (i=0; i<Get_nElems(); i++) */
      
  return val;
}
//_____________________________________________________________________________
double NPoly::Eval(const vecd &X) const 
{
  //same as above, but use the coefficients 'hard-coded' to each element
  
  if ((int)X.size() != Get_nDoF()) { 
    Error("Eval(vecd)", "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), (int)Get_nDoF());
    return -1e30;
  }
      
  //now, actually evaluate
  double val=0.; 

  const int max_pow = Get_maxPower(); 
  double X_pows[Get_nDoF()][max_pow + 1]; 
  for (int d=0; d<Get_nDoF(); d++) {
    for (int p=0; p<=max_pow; p++) X_pows[d][p] = pow(X[d], p); 
  }

  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff.
  for (const auto &elem : fElems) {
    
    double elem_val=1.;
    
    for (int d=0; d<Get_nDoF(); d++) {
      if (elem.powers[d]>0) elem_val *= X_pows[d][elem.powers[d]];
    }
    //check to see if this is the one element for which all pows. are zero
    val += elem_val * elem.coeff; 
	
  }//for (const auto &elem : fElems)
      
  return val;
}
//_____________________________________________________________________________
vecd NPoly::Gradient(const vecd &coeff, const vecd &X) const
{
  if (coeff.size() != Get_nElems()) { 
    Error("NPoly::Eval()",
	  "Wrong number of coeffs. given; expected %i, recieved %i.",
	  Get_nElems(), (int)coeff.size() );
    return {};
  }
  
  if ((int)X.size() != fnDoF) { 
    Error("NPoly::Eval_noCoeff",
	  "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), fnDoF);
    return {};
  }

  vecd grad(fnDoF,0.); 
  
  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff. 
  for (unsigned int i=0; i<Get_nElems(); i++) { auto elem = fElems[i];
    
    vecd grad_elem(fnDoF,1.); 

    for (int g=0; g<fnDoF; g++) { 

      int pow_deriv = elem.powers[g];
      if (pow_deriv<1) {
        grad_elem[g] = 0.; continue; 
      } else { 
        grad_elem[g] *= ((double)pow_deriv)*pow(X[g], pow_deriv-1);
      } 

      //now, deal with the elems we are NOT taking the deriv of 
      for (int d=0; d<fnDoF; d++) { if (d==g) continue; 
        if (elem.powers.at(d)>0) { 
          grad_elem[g] *= pow( X.at(d), elem.powers.at(d) );
        }
      }
    }

    
    grad += coeff.at(i) * grad_elem; 
  }
  
  //for (i=0; i<Get_nElems(); i++) 
  return grad;
}
//_____________________________________________________________________________
vecd NPoly::Gradient(const vecd &X) const 
{
  //same as above, but use the coefficients 'hard-coded' to each element
  
  if ((int)X.size() != Get_nDoF()) { 
    Error("Gradient()", "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), (int)Get_nDoF());
    return {};
  }
  
  //now, actually evaluate
  RVec<double> ret(Get_nDoF(), 0.); 

  const int max_pow = Get_maxPower(); 
  double X_pows[Get_nDoF()][max_pow + 1]; 
  for (int d=0; d<Get_nDoF(); d++) {
    for (int p=0; p<=max_pow; p++) X_pows[d][p] = pow(X[d], p); 
  }

  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff.
  for (const auto &elem : fElems) {
    
    double elem_val=1.;
    
    for (int d=0; d<Get_nDoF(); d++) {
      if (elem.powers[d]>0) elem_val *= X_pows[d][elem.powers[d]];
    }
    //multiply by the coefficient of this element
    elem_val *= elem.coeff; 

    //use the power rule to take the derivative w/r/t each input variable
    for (int d=0; d<Get_nDoF(); d++) {
      if (elem.powers[d] > 0) 
        ret[d] += elem_val * ((double)elem.powers[d]) / X[d]; 
    }
  }//for (const auto &elem : fElems)
      
  return ret;
}
//_____________________________________________________________________________
RMatrix NPoly::Hessian(const vecd &X) const 
{
  //same as above, but use the coefficients 'hard-coded' to each element
  
  if ((int)X.size() != Get_nDoF()) { 
    Error("Hessian()", "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), (int)Get_nDoF());
    return {};
  }
  
  //now, actually evaluate
  RMatrix ret(Get_nDoF(), Get_nDoF(), 0.); 

  const int max_pow = Get_maxPower(); 
  double X_pows[Get_nDoF()][max_pow + 1]; 
  for (int d=0; d<Get_nDoF(); d++) {
    for (int p=0; p<=max_pow; p++) X_pows[d][p] = pow(X[d], p); 
  }

  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff.
  for (const auto &elem : fElems) {
    
    double elem_val=1.;
    
    for (int d=0; d<Get_nDoF(); d++) {
      if (elem.powers[d]>0) elem_val *= X_pows[d][elem.powers[d]];
    }
    //multiply by the coefficient of this element
    elem_val *= elem.coeff; 

    //use the power rule to take the derivative w/r/t each input variable
    for (int di=0; di<Get_nDoF(); di++) {
      for (int dj=di; dj<Get_nDoF(); dj++) {
        
        double val(elem_val / (X[di] * X[dj])); 

        if (di==dj) { val *= (double)(elem.powers[di] * (elem.powers[di]-1)); }
        else        { val *= (double)(elem.powers[di] * elem.powers[dj]);     }
        
        ret.get(di,dj) += val; 
      }
    }
  }//for (const auto &elem : fElems)

  //because the hessian is symmetric, we only computed the elements once. now we must fill in their 'mirror' elements
  for (int di=1; di<Get_nDoF(); di++) 
    for (int dj=0; dj<di; dj++) 
      ret.get(di,dj) = ret.get(dj,di); 

      
  return ret;
}//_____________________________________________________________________________
RVec<int> NPoly::Get_elemPowers(unsigned int i) const
{
  if (i >= Get_nElems()) {
    Warning("NPoly::Get_elemPowers",
	    "Invalid element-index requested; %i, max is %i",i,Get_nElems()-1);
    return {};
  }
  return fElems.at(i).powers;
}
//_____________________________________________________________________________
double NPoly::Get_elemCoeff(unsigned int i) const
{
  if (i >= Get_nElems()) {
    Warning("Get_elemCoeff",
	    "Invalid element-index requested; %i, max is %i",i,Get_nElems()-1);
    return {};
  }
  return fElems.at(i).coeff;
}
//_____________________________________________________________________________
const NPoly::NPolyElem* NPoly::Get_elem(unsigned int i) const 
{
  if (i >= Get_nElems()) {
    Warning("Get_elemCoeff",
	    "Invalid element-index requested; %i, max is %i",i,Get_nElems()-1);
    return nullptr;
  }
  return &fElems[i];
}
//_____________________________________________________________________________
void NPoly::Add_element(const RVec<int> &pows, double coefficient)
{ 
  //just calls the method below. 
  this->Add_element({.powers=pows, .coeff=coefficient}); 
}
//_____________________________________________________________________________
void NPoly::Add_element(const NPoly::NPolyElem& elem)
{
  if (elem.powers.size()!=Get_nDoF()) {
    Error("Add_element", "Num. of powers (%i) incorrect size for this Poly's DoF (%i)",
	  (int)elem.powers.size(), Get_nDoF());
    return;
  }
  
  //for the qucikened version of Eval(), we must know what the maximum exponential power is for any input var. 
  //inscrutable segfaults will result if not!!
  for (const int& pow : elem.powers) if (pow > f_maxPower) f_maxPower=pow; 

  fElems.push_back(elem); 
}
//_____________________________________________________________________________
void NPoly::AutoConstructPoly(const int max_power, const int nDoF)
{
  RVec<int> powers(max_power,0); 
  
  while (1) {
    //create a new poly element with our current 'powers'
    //auto elem = NPolyElem;
    RVec<int> pows; 
    
    for (int x=0; x<nDoF; x++) {
      int myPow=0;
      for (int xPow : powers) if (xPow==x) myPow++;
      
      pows.push_back(myPow);
    }
    
    /*cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
      for (int p=0; p<max_power; p++) { 
      for (int d=0; d<nDoF+1; d++) cout << (powers.at(p)==d?1:0) << " ";
      cout << endl; 
      }*/ 

    /*cout << "elem : "; 
      for (int d=0; d<nDoF; d++) cout << elem->powers.at(d) << ","; 
      cout << endl; */ 

    this->Add_element(pows, 1.);
    //fElems.push_back({.powers=pows, .coeff=1.});
    
    //now, increase to the next set of powers
    for (int it1=powers.size()-1; it1>=0; it1--) { 
      //go back up and reset higher powers
            
      if (powers.at(it1)<nDoF) {
        powers.at(it1) += 1; 
        for (unsigned int it2=it1+1; it2<powers.size(); it2++) {
          powers.at(it2) = powers.at(it1);
        }
        break; 
      }
      
      if (it1==0) {
        //cout << "polynomial size = " << elems.size() << endl; 
        return; 
      }
            
    }//for (int it1=powers.size()-1 ... [power-escalating group]        
  }//while (1)
}
//_____________________________________________________________________________
void NPoly::Print() const 
{
  printf("Polynomial DoF: %i, all elements (%i)\n", Get_nDoF(), Get_nElems()); 

  for (int i=0; i<Get_nElems(); i++) {
    
    const auto elem = Get_elem(i); 

    printf(" -- {"); 
    
    for (int j=0; j<elem->powers.size(); j++) printf(" %2i ", elem->powers.at(j)); 
    
    printf("}  -- %+.8e\n", elem->coeff); 
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
//these following methods will be used for symbolic computation; 
// for example, symbolically computing the result of feeding the output of one NPolyArray into the input of another. 
// I am experimenting with constructing optics models which are composed of several NPolyArray's 'chained together' 
// in which the output of each NPolyArray in the chain is fed into the input of the next. 
//_____________________________________________________________________________
NPoly NPoly::operator*(const NPoly& rhs) const
{
  if (rhs.Get_nDoF() != Get_nDoF()) {
    Error("operator*", "LHS NPoly DoF (%i) does not match RHS NPoly DoF(%i)",
	  rhs.Get_nDoF(), Get_nDoF());
    return NPoly(0);
  }

  NPoly prod(Get_nDoF()); 

  vector<NPoly::NPolyElem> new_elems; new_elems.reserve(Get_nElems() * rhs.Get_nElems());  


  //now, just do the 'foil' of multiplying each element together. 
  for (int i=0; i<Get_nElems(); i++) {
    for (int j=0; j<rhs.Get_nElems(); j++) {

      //recalling how to do this from grade school...

      //exponents add together
      RVec<int> pows = Get_elemPowers(i) + rhs.Get_elemPowers(j);
      
      //coefficents multiply 
      double coeff = Get_elemCoeff(i) * rhs.Get_elemCoeff(j);
      
      //now, add this new element to the set of all elements which will constitute our new polynomial. 
      new_elems.push_back({.powers=pows, .coeff=coeff}); 
    }
  }

  //now, to make our lives easier down the road, we will see if any of the new elements we just constructed above
  // can be added together. this is done by checking for elements which have idential 'exponent signatures', i.e. all the same exponents.  
  vector<NPoly::NPolyElem> new_elems_unique{}; 

  for (int i=0; i<new_elems.size(); i++) {

    auto elem = new_elems.at(i); 

    //check to see if an element just like this one exists in the vector already. 
    auto it = find( new_elems_unique.begin(), new_elems_unique.end(), elem ); 
    
    if (it==new_elems_unique.end()) {//an element with the same set of exponents does NOT yet exist in the 'new_elems_unique' vec

      new_elems_unique.push_back(elem); 
      
    } else { //if we've gotten here, it means that there is an element in the 'new_elems_unique' vector with the same exponents. 

      //so we can 'combine like terms' here. 
      //think of adding: 
      //    P1(x) = x + y^2;
      //    P2(x) = 3 y^2 
      //    (P1 + P2)(x) = x + 4 y^2. 
      
      //this iterator now points to the elem in 'new_elems_unique' which is a match 
      it->coeff += elem.coeff;
    }
  }

  //now that we have combined like terms, we can add all our elements to our vec. 
  for (const auto& elem : new_elems_unique) prod.Add_element(elem); 

  return prod;  
}
//_____________________________________________________________________________
NPoly NPoly::Pow(const NPoly& pol, const int pow) 
{
  if (pow < 0) {
    fprintf(stderr, "Error in <NPoly::Pow(NPoly&, int)>: Power %i passed as exponent, only non-negative integers are supported.", pow); 
    return NPoly(0); 
  }

  NPoly ret(pol.Get_nDoF()); 

  //we add the element to this polynomal '1.', so that we can just multiply it by 'pol' as many times as we need to 
  // raise it to the power requested. 
  ret.Add_element(RVec<int>(pol.Get_nDoF(), 0.), 1.);  

  for (int i=0; i<pow; i++) ret = ret * pol; 

  return ret; 
}
//_____________________________________________________________________________
NPoly NPoly::operator+(const NPoly& rhs) const
{
  if (rhs.Get_nDoF() != Get_nDoF()) {
    Error("operator+", "LHS NPoly DoF (%i) does not match RHS NPoly DoF(%i)",
	  Get_nDoF(), rhs.Get_nDoF());
    return NPoly(0);
  }

  NPoly sum(Get_nDoF()); 

  //add all elements to this vector. 
  vector<NPoly::NPolyElem> new_elems; new_elems.reserve(Get_nElems() + rhs.Get_nElems());  

  for (int i=0; i<Get_nElems(); i++) { 
    new_elems.push_back(*(Get_elem(i))); 
  }

  for (int i=0; i<rhs.Get_nElems(); i++) { 
    new_elems.push_back(*(rhs.Get_elem(i))); 
  }

  //now combine like terms to eliminate redundant elements. 
  vector<NPoly::NPolyElem> new_elems_unique{}; 

  for (int i=0; i<new_elems.size(); i++) {

    auto elem = new_elems.at(i); 

    //check to see if an element just like this one exists in the vector already. 
    auto it = find( new_elems_unique.begin(), new_elems_unique.end(), elem ); 
    
    if (it==new_elems_unique.end()) {//an element with the same set of exponents does NOT yet exist in the 'new_elems_unique' vec

      new_elems_unique.push_back(elem); 
      
    } else { //if we've gotten here, it means that there is an element in the 'new_elems_unique' vector with the same exponents. 

      //so we can 'combine like terms' here. 
      //think of adding: 
      //    P1(x) = x + y^2;
      //    P2(x) = 3 y^2 
      //    (P1 + P2)(x) = x + 4 y^2. 
      
      //this iterator now points to the elem in 'new_elems_unique' which is a match 
      it->coeff += elem.coeff;
    }
  }

  //now that we have combined like terms, we can add all our elements to our vec. 
  for (const auto& elem : new_elems_unique) sum.Add_element(elem); 

  return sum; 

}
//_____________________________________________________________________________

ClassImp(NPoly)

