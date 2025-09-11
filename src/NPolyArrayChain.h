#ifndef NPolyArrayChain_h_
#define NPolyArrayChain_h_

//
//  Chain of NPolyArray objects, with the output of one feeding the next. 
//  
//
//
//

#include <vector>
#include "NPolyArray.h"
#include <utility> 
#include <ROOT/RVec.hxx> 
#include <RMatrix.h> 


class NPolyArrayChain {

public: 
    //this is the allowed set of types that a 'link' in the polynomial chain can be. 
    enum EArrayType { 
        kStatic =0,     //'static', i.e., the elements of this polynomial will not be modified.  
        kMutable =1,    //'mutable', i.e., the elements of this polynomial *can* be modified. 
        kBuffer =2,     //this is a 'buffer' polynomial, which will be 'folded into' a polynomial on either side of it. 
    }; 

    NPolyArrayChain() : arrays{} {}; 
    ~NPolyArrayChain() {}; 
    
    //add an array on top of the last output array. 
    void AppendArray(NPolyArray arr, bool is_mutable=NPolyArrayChain::kStatic); 

    //add a 'buffer' polynomial. 
    void InsertBufferArray( const int DoF, const int order=1 ); 


    //evaluate all arrays, starting with arrays[0], feeding the input of each one into the next. 
    ROOT::RVec<double> Eval(const ROOT::RVec<double>& X) const;
    RMatrix Jacobian(const ROOT::RVec<double>& X) const; 

    int Get_DoF_in()  const { return arrays.front().first.Get_DoF_in(); }
    int Get_DoF_out() const { return arrays.back().first.Get_DoF_out(); }

    size_t N_arrays() const { return arrays.size(); }

    int Get_nElems_mutable() const; 

    //add an update to each of the weights. 
    NPolyArrayChain& operator+=(const ROOT::RVec<double>& dW); 

    //Returns an NPolyArray of the given DoF and oder, with only the 'linear' elements nonzero (in other words, this poly array is the identity matrix). 
    static NPolyArray CreateIdentityNPolyArray(const int DoF, const int order); 

    //the first member of the pair is the array, and the second member of the array is 'true' if the elements of this array 
    // can be modified, and 'false' if they cannot. 

    std::vector<std::pair<NPolyArray, NPolyArrayChain::EArrayType>> arrays; 

private: 
    
    ClassDef(NPolyArrayChain, 1); 
}; 


#endif 