
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

using namespace std; 

//this will be a temporary funct., that I may absorb into NPoly.cxx. 
int parse_poly_from_file(const char* path_dbfile, NPoly *poly) 
{
    return 0; //noop
}

int test_forward_model( const char* path_infile="",
                        const char* path_dbfile="data/csv/db_test.dat",  
                        const char* tree_name="track_data" ) 
{
    const char* const here = "test_forward_model"; 

    auto infile = new TFile(path_infile, "READ");

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(here, "libApexOptics could not be loaded."); 
        return 1; 
    }
    
    //check if the infile could be opened
    if (!infile || infile->IsZombie()) {
        Error(here, "root file '%s' could not be opened.", path_infile); 
        return 1; 
    }

    //check if we can find the proper TTree
    TTree* tree = (TTree*)infile->Get(tree_name); 
    if (!tree) {
        Error(here, "could not find TTree '%s'", tree_name); 
        return 1; 
    }
       
    delete tree; 
    infile->Close(); 
    delete infile; 


    //check if we can open data file
    fstream dbfile(path_dbfile, ios::in); 
    
    if (!dbfile.is_open()) {
        Error(here, "could not open db file '%s'", path_dbfile); 
        return 1; 
    }


    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    return 0; 
}