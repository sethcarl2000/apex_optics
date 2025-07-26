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
    NPolyArray(int _DoF_out, int _DoF_in);
    NPolyArray(const std::vector<NPoly>& _polys);

    ~NPolyArray() {}; 

    inline int Get_DoF_in()  const { return fDoF_in;  }
    inline int Get_DoF_out() const { return fDoF_out; }

    ROOT::RVec<double> Eval(const ROOT::RVec<double> &X) const; 

    //Each row is the gradient of one of the constituent polynomials
    RMatrix Jacobian(const ROOT::RVec<double> &X) const;
     
    //Each element is the hessian of one of the polynomials. 
    ROOT::RVec<RMatrix> HessianTensor(const ROOT::RVec<double>& X) const; 

    //get the underlying polynomial which is responsible for a coordinate of this array. 
    inline NPoly* Get_poly(int i); 

    //basically a const_cast-ed version of the above. 
    inline const NPoly* Get_poly_const(int i) const { return Get_poly(i); }

    //this symbolically computes the action of feeding the output of one NPolyArray into the input of another ('nesting'). 
    // This way, if an optics model is composed of a 'chain' of NPolyArray's, in which the output of one is fed into the input of the next, 
    // we will be able to symbolically compute that action once, and not have to numerically compute each polynomial in the chain separately.
    // this has the potential to reduce a large number of redundant calculations, as well as making the logic needed to compute Jacobians/Hessians
    // much simpler. 
    static NPolyArray Nest(const NPolyArray& output, const NPolyArray &input); 


private: 

    int fDoF_out, fDoF_in; 

    std::vector<NPoly> fPolys; 

    ClassDef(NPolyArray,1);
};


#endif