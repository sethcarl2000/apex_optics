
#include "NPolyArray.h"
#include "RMatrix.h"
#include <algorithm>
#include <stdio.h>
#include <cstdio> 


using namespace std; 
using namespace ROOT::VecOps; 

//______________________________________________________________________________________________
NPolyArray::NPolyArray(int _DoF_out, int _DoF_in) 
    : fDoF_out(_DoF_out), fDoF_in(_DoF_in)
{
    //create empty polynomials with the proper dimensions
    for (int i=0; i<Get_DoF_out(); i++) 
        fPolys.push_back(NPoly(Get_DoF_in()));
}

//______________________________________________________________________________________________
NPolyArray::NPolyArray(const vector<NPoly>& _polys) 
{
    const char* const here = "NPolyArray(vector<NPoly>)";
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
NPolyArray::NPolyArray(const vector<NPoly*>& _polys) 
{
    const char* const here = "NPolyArray(vector<NPoly>)";
    //use the first polynomial to check the number of input DoF
    fDoF_out = _polys.size(); 
    fDoF_in  = 0; 
    
    if (fDoF_out < 1) { 
        Error(here, "Input vector is empty. Cannot construct"); 
        return; 
    }

    //check to make sure that each poly has the same input DoF
    fDoF_in = _polys.at(0)->Get_nDoF(); 
    for (const NPoly* poly : _polys) {
        if (fDoF_in != poly->Get_nDoF()) {
            Error(here, "DoF of input polys does not match. Cannot construct"); 
            fDoF_in  =0; 
            fDoF_out =0; 
            return; 
        } 
    }

    //if we got here, then all polys match one another's DoF. We may now construct our vector of NPolys 
    fPolys.reserve(_polys.size()); 
    for (const NPoly* poly : _polys) fPolys.push_back(*poly); 
}

//______________________________________________________________________________________________
RVec<double> NPolyArray::Eval(const RVec<double>& X) const 
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
const NPoly *NPolyArray::Get_poly(int i) const 
{
    if (i<0 || i>=Get_DoF_out()) {
        Error("Get_poly(int)", "Acess to array of polynomals out-of-range; asked for index %i, range is [0,%i].", 
                i, Get_DoF_out()-1);
        return nullptr; 
    }
    return &(fPolys.at(i)); 
}
//______________________________________________________________________________________________
RMatrix NPolyArray::Jacobian(const RVec<double> &X) const 
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
RVec<RMatrix> NPolyArray::HessianTensor(const RVec<double> &X) const 
{
    if ((int)X.size() != Get_DoF_in()) {
        Error("HessianTensor()", "Input vector wrong size; got %i, expected %i", (int)X.size(), Get_DoF_in());
        return {}; 
    }
    
    vector<RMatrix> ret; ret.reserve(Get_DoF_out()); 

    for (int i=0; i<Get_DoF_out(); i++) ret.push_back(Get_poly(i)->Hessian(X)); 
    
    return ret; 
}
//______________________________________________________________________________________________
NPolyArray NPolyArray::Nest(const NPolyArray& output, const NPolyArray& input) 
{   
    //the number of inputs to the 'output' array must match the number of outputs to the 'input' array
    if (output.Get_DoF_in() != input.Get_DoF_out()) {
        fprintf(stderr, "Error in <NPolyArray::Nest>: The input-DoF of the output array (%i) does not match the output-DoF of the input array (%i)\n", 
        output.Get_DoF_in(), input.Get_DoF_out()); 
        return NPolyArray(0,0); 
    }

    //make a copy of the inputs' polys
    vector<NPoly> input_polys{}; 
    for (int i=0; i<input.Get_DoF_out(); i++) input_polys.push_back( *(input.Get_poly(i)) ); 

    vector<NPoly> output_polys;

    for (int i=0; i<output.Get_DoF_out(); i++) {
        
        NPoly out_nested(output.Get_DoF_in()); 

        const NPoly *out_pol = output.Get_poly(i); 

        for (int e=0; e<out_pol->Get_nElems(); e++) {

            const NPoly::NPolyElem* out_pol_elem = out_pol->Get_elem(e); 
            
            //make a new NPoly which only has one element = 1. 
            NPoly out_nested_elem(output.Get_DoF_in()); 
            out_nested_elem.Add_element(RVec<int>(output.Get_DoF_in(), 0.), 1.); 

            //multiply all the input polynomials together, raising each to the power proscribed by the output polynomials. 
            for (int j=0; j<output.Get_DoF_in(); j++) { 
                out_nested_elem = out_nested_elem * NPoly::Pow( input_polys.at(j), out_pol_elem->powers.at(j) ); 
            }
            
            out_nested_elem *= out_pol_elem->coeff; 

            //add this 'element' to the output polynomial 
            out_nested = out_nested + out_nested_elem; 
        }

        output_polys.push_back( out_nested ); 
    }
    
    return NPolyArray(output_polys); 
}
//______________________________________________________________________________________________
//______________________________________________________________________________________________

ClassImp(NPolyArray); 