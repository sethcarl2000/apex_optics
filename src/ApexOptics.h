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
#include <map>

namespace ApexOptics {

    //given an array of input branches, and a 'target' output branch, will attempt to 
    // fit a polynomial (of order poly_order) from the inputs to the output. 
    std::map<std::string, NPoly*> Create_NPoly_fit( ROOT::RDF::RNode df, 
                                                    const int poly_order, 
                                                    const std::vector<std::string> &inputs, 
                                                    const std::vector<std::string> &outputs ); 

    //Given an std::map<string,NPoly*>, where the 'string' key is the name of the polynomial to be written to the 
    //file, this fcn will create (and truncate!) a '.dat' file which contains all information needed to 
    //reconstruct an NPoly's elements. 
    int Create_dbfile_from_polymap(bool is_RHRS, std::string path_outfile, std::map<std::string, NPoly*> polymap); 

    //Given a path to a db-file to use, the polynomial name, and a ptr to a polynomial to use, this will fill 
    // all NPoly elements found in the dbfile into the polynomial. returns the number of elements found. 
    // return code < 0 means that something has failed in the opening of the dbfile. 
    int Parse_NPoly_from_file(const char* path_dbfile, const char* poly_name, NPoly *poly);
};

#endif