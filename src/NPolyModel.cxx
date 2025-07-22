
#include "NPolyModel.h"
#include "RMatrix.h"
#include <algorithm>
#include <stdio.h>


using namespace std; 
using namespace ROOT::VecOps; 

//______________________________________________________________________________________________
NPolyModel::NPolyModel(int _DoF_in, int _DoF_out) 
    : fDoF_in(_DoF_in), fDoF_out(_DoF_out)
{
    //create empty polynomials with the proper dimensions
    for (int i=0; i<Get_DoF_out(); i++) 
        fPolys.push_back(NPoly(Get_DoF_in()));
}
//______________________________________________________________________________________________
NPolyModel::NPolyModel(const vector<NPoly>& _polys) 
{
    const char* const here = "NPolyModel(vector<NPoly>)";
    //use the first polynomial to check the number of input DoF
    fDoF_out = _polys.size(); 
    fDoF_in  = 0; 
    
    if (fDoF_out < 1) { 
        Error(here, "Input vector is empty. Cannot construct"); 
        return; 
    }

    fDoF_in = _polys.at(0).Get_nDoF(); 
    for (const NPoly& poly : _polys) {
        if (fDoF_in != poly.Get_nDoF()) {
            Error(here, "DoF of input polys does not match. Cannot construct"); 
            fDoF_in  =0; 
            fDoF_out =0; 
            return; 
        } 
    }

    //if we got here, then all polys match one another. 
    fPolys = _polys; 
}
//______________________________________________________________________________________________
RVec<double> NPolyModel::Eval(const RVec<double>& X) const 
{
    const char* const here = "Eval(const RVec<double>& X)"; 

    if ((int)X.size() != Get_DoF_in()) {
        Error(here, "Input vec wrong size; %i, expected %i.", (int)X.size(), Get_DoF_in()); 
        return {}; 
    }

    RVec<double> ret(Get_DoF_out(), 0.); 

    for (int i=0; i<Get_DoF_out(); i++) ret[i] = fPolys[i].Eval(X); 

    return ret; 
}
//______________________________________________________________________________________________
NPoly *NPolyModel::Get_poly(int i)
{
    if (i<0 || i>=Get_DoF_out()) {
        Error("Get_poly(int)", "Acess to array of polynomals out-of-range; asked for index %i, range is [0,%i].", 
                i, Get_DoF_out()-1);
        return nullptr; 
    }
    return &(fPolys.at(i)); 
}
//______________________________________________________________________________________________
RMatrix NPolyModel::Jacobian(const RVec<double> &X) const 
{
    if ((int)X.size() != Get_DoF_in()) {
        Error("Jacobian()", "Input vector wrong size; got %i, expected %i", (int)X.size(), Get_DoF_in());
        return RMatrix(0,0); 
    }
    
    vector<double> J_vec; 
    J_vec.reserve(Get_DoF_in()*Get_DoF_out()); 

    //compute the gradient of each poly, then add the results into a single vector
    for (int i=0; i<Get_DoF_out(); i++) {
        
        RVec<double> grad{fPolys[i].Gradient(X)}; 

        //this may be a bit faster than going thru and copying each element manually
        copy(grad.begin(), grad.end(), back_inserter(J_vec)); 
    }

    //now, create the matrix
    return RMatrix(Get_DoF_out(), Get_DoF_in(), J_vec); 
}
//______________________________________________________________________________________________
//______________________________________________________________________________________________
//______________________________________________________________________________________________
//______________________________________________________________________________________________

ClassImp(NPolyModel); 