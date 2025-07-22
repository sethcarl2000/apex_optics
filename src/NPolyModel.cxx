
#include "NPolyModel.h"
#include "RMatrix.h"

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
    const char* const here =  "NPolyModel(vector<NPoly>)";
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
//______________________________________________________________________________________________
//______________________________________________________________________________________________
//______________________________________________________________________________________________
//______________________________________________________________________________________________
//______________________________________________________________________________________________
//______________________________________________________________________________________________
