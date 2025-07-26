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
    NPolyArray(int _DoF_in, int _DoF_out);
    NPolyArray(const std::vector<NPoly>& _polys);

    ~NPolyArray() {}; 

    inline int Get_DoF_in()  const { return fDoF_in;  }
    inline int Get_DoF_out() const { return fDoF_out; }

    ROOT::RVec<double> Eval(const ROOT::RVec<double> &X) const; 

    RMatrix Jacobian(const ROOT::RVec<double> &X) const; 

    NPoly* Get_poly(int i); 

private: 

    int fDoF_in, fDoF_out; 

    std::vector<NPoly> fPolys; 

    ClassDef(NPolyArray,1);
};


#endif