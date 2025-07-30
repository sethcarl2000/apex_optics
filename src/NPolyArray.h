#ifndef NPolyArray_h_
#define NPolyArray_h_

//////////////////////////////////////////////////////////////////////////
//
//  NPolyArray
// 
//  This is meant to be an Optics model which maps from a N-space to an M-space
//  each coordinate in M-space represented by a disticnt NPoly. 
//
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TObject.h"
#include "NPoly.h"
#include "RMatrix.h"
#include <ROOT/RVec.hxx> 
#include <vector>

class NPolyArray : public TObject {

public: 
    
    //default constructor (NPolys empty-initialized with correct DoF)
    NPolyArray(int _DoF_out, int _DoF_in);
    
    //constructor with vector of NPoly's (by value)
    NPolyArray(const std::vector<NPoly>& _polys);

    //constructor with vector of NPoly's (by ptr)
    NPolyArray(const std::vector<NPoly*>& _polys); 

    ~NPolyArray() {}; 

    inline int Get_DoF_in()  const { return fDoF_in;  }
    inline int Get_DoF_out() const { return fDoF_out; }

    ROOT::RVec<double> Eval(const ROOT::RVec<double> &X) const; 

    //Each row is the gradient of one of the constituent polynomials
    RMatrix Jacobian(const ROOT::RVec<double> &X) const;
     
    //Each element is the hessian of one of the polynomials. 
    ROOT::RVec<RMatrix> HessianTensor(const ROOT::RVec<double>& X) const; 

    //get the underlying polynomial which is responsible for a coordinate of this array. 
    const NPoly* Get_poly(int i) const; 


    //this symbolically computes the action of feeding the output of one NPolyArray into the input of another ('nesting'). 
    // This way, if an optics model is composed of a 'chain' of NPolyArray's, in which the output of one is fed into the input of the next, 
    // we will be able to symbolically compute that action once, and not have to numerically compute each polynomial in the chain separately.
    // this has the potential to reduce a large number of redundant calculations, as well as making the logic needed to compute Jacobians/Hessians
    // much simpler. 
    static NPolyArray Nest(const NPolyArray& output, const NPolyArray &input); 


    //what we're effectivley doing here is using newton's method for iterating towrads the root of a nonlinear system.
    // the system we're trying to solve is the following least-square problem: 
    //  - the 'chi-square' in this case is the square error between the model's value for Xfp, and the actual value.
    //    this is given by the d_Xfp vector above. 
    //  - Therefore, the 'F' vector is our evaluation of the gradient of this function, which will be zero at the 
    //    minimum error value (if it isn't a local, false minima.)
    //  - to use newton's method, we need to compute the Jacobian of our 'F' funciton. this is what 'J' will be. 
    //
    int Iterate_to_root(ROOT::RVec<double>& X, const ROOT::RVec<double>& Z, int n_iterations=1) const;


private: 

    int fDoF_out, fDoF_in; 

    std::vector<NPoly> fPolys; 

    ClassDef(NPolyArray,1);
};


#endif