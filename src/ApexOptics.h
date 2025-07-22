#ifndef ApexOptics_h_
#define ApexOptics_h_

//////////////////////////////////////////////////////////////////////////
//
// ApexOptics
//
// This is a namespace defintion to give us some helper functions which 
// can be used in conjuction with the ApexOptics library. 
//
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include <ROOT/RDataFrame.hxx>  
#include "NPoly.h"
#include <vector>
#include <string>
#include <memory>

namespace ApexOptics {

    //given an array of input branches, and a 'target' output branch, will attempt to 
    // fit a polynomial (of order poly_order) from the inputs to the output. 
    NPoly* Create_NPoly_fit(ROOT::RDF::RNode df, 
                            const int poly_order, 
                            const std::vector<std::string> &inputs, 
                            const char* output); 

};

#endif