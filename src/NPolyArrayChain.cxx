#include "NPolyArrayChain.h"
#include "NPoly.h"
#include <sstream> 
#include <stdexcept> 

using namespace std; 
using namespace ROOT::VecOps; 

//____________________________________________________________________________________________________________________________
void NPolyArrayChain::AppendArray(NPolyArray arr, bool is_mutable) { 

    //check if this array matches with the others already added
    if (!arrays.empty()) {
        if (arr.Get_DoF_in() != arrays.back().first.Get_DoF_out()) {
            ostringstream oss; 
            oss << "in <NPolyArrayChain::AppendArray>: input DoF of new array (" << arr.Get_DoF_in() << ")"
                    " does not match output DoF of array it will attach to (" << arrays.back().first.Get_DoF_out() << ")"; 
            throw invalid_argument(oss.str());  
            return; 
        }
    }
    arrays.emplace_back( arr, (is_mutable ? kMutable : kStatic) ); 
}
//____________________________________________________________________________________________________________________________
void NPolyArrayChain::InsertBufferArray( const int DoF, const int order ) {

    //check if this array matches with the others already added
    if (!arrays.empty()) {
        if (DoF != arrays.back().first.Get_DoF_out()) {
            ostringstream oss; 
            oss << "in <NPolyArrayChain::InsertBufferArray>: input DoF of new array (" << DoF << ")"
                    " does not match output DoF of array it will attach to (" << arrays.back().first.Get_DoF_out() << ")"; 
            throw invalid_argument(oss.str());  
            return; 
        }
    }
    arrays.emplace_back( NPolyArrayChain::CreateIdentityNPolyArray(DoF, order), kBuffer ); 
}
//____________________________________________________________________________________________________________________________
//evaluate all arrays, starting with arrays[0], feeding the input of each one into the next. 
RVec<double> NPolyArrayChain::Eval(const RVec<double>& X) const { 
    RVec<double> out{X}; 
    for (const auto& array : arrays) out = std::move( array.first.Eval(out) );
    return out;  
}
//____________________________________________________________________________________________________________________________
int NPolyArrayChain::Get_nElems_mutable() const {
    //get the total number of 'mutable' parameters in the polyarraychain
    int n_elems=0; 
    for (const auto& array : arrays) { if (array.second != kStatic) n_elems += array.first.Get_nElems(); } 
    return n_elems; 
}
//____________________________________________________________________________________________________________________________
NPolyArrayChain& NPolyArrayChain::operator+=(const RVec<double>& dW) 
{
    //check if the input vector is the right size 
    if (dW.size() != (size_t)Get_nElems_mutable()) {
        ostringstream oss; 
        oss << "in <NPolyArrayChain::operator+=>: number of mutable weights (" << Get_nElems_mutable() << ")"
                " does not match size of input RVec<double> (" << dW.size() << ")"; 
        throw std::invalid_argument(oss.str()); 
        return *this; 
    }

    //now, add each element of 'dW' to our existing 'mutable' weights
    //save the current weights as the best yet found
    int i_elem=0;   
    for (int i=N_arrays()-1; i>=0; i--) { 
        
        //check if this array is marked as mutable
        if (arrays[i].second == kStatic) continue; 

        auto& parr = arrays[i].first; 

        for (int j=0; j<parr.Get_DoF_out(); j++) { 
            
            auto poly = parr.Get_poly(j); 
            for (int k=0; k<poly->Get_nElems(); k++) poly->Get_elem(k)->coeff += dW.at(i_elem++);
        } 
    }//for (int i=sandwich->N_arrays()-1; i>=0; i--)

    return *this; 
}
//____________________________________________________________________________________________________________________________
NPolyArray NPolyArrayChain::CreateIdentityNPolyArray(const int DoF, const int order)
{   
    vector<NPoly> polys; 
    for (int i=0; i<DoF; i++) {

        NPoly poly(DoF, order); 
        //set all elements to zero 
        for (int j=0; j<poly.Get_nElems(); j++) poly.Get_elem(j)->coeff = 0.; 

        //Now, set all 'linear' elments to 1.
        RVec<int> powers(DoF, 0); 
        powers[i] = 1; 

        poly.Find_element(powers)->coeff = 1.; 

        polys.push_back(poly); 
    }

    return NPolyArray(polys); 
};
//____________________________________________________________________________________________________________________________

ClassImp(NPolyArrayChain); 