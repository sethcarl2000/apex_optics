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
  
  for (const NPolyElem &elem : fElems) {

    double elem_val=1.;

    for (int d=0; d<Get_nDoF(); d++) {
      if (elem.powers[d] > 0) elem_val *= pow( X[d], elem.powers[d] );
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
    Error("Eval(vecd)", "Size of input vector (%i) does not match poly nDoF (%i)",
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
//_____________________________________________________________________________
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
NPoly::NPolyElem* NPoly::Get_elem(unsigned int i) 
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
  if (pows.size()!=fnDoF) {
    Error("Add_element", "Num. of powers (%i) incorrect size for this Poly's DoF (%i)",
	  (int)pows.size(), (int)fnDoF);
    return;
  }
  //record the maximum exponent power in this polynomial 
  for (const int& pow : pows) if (pow > f_maxPower) f_maxPower=pow; 

  fElems.push_back({.powers=pows, .coeff=coefficient}); 
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
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________

ClassImp(NPoly)

