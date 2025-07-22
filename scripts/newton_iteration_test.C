#include "TROOT.h"

using namespace std; 


//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
int fitpoints_mc_sv_fp(  bool is_RHRS=false,
                         const int poly_order=2,
                         const char* path_infile="",
                         const char* path_dbfile="data/csv/db_prod_mc_sv_fp_L_2ord.dat",  
                         const char* tree_name="tracks_fp" ) 
{
    const char* const here = "fitpoints_mc_sv_fp"; 

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

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 


    const double hrs_momentum = 1104.0; 

    vector<string> branches_fp = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    const int DoF_sv = 5; 

    vector<NPoly> poly_vec; 

    for (const string& str : branches_fp) {

        NPoly poly(DoF_sv); 

        ApexOptics::Parse_NPoly_from_file(path_dbfile, str.data(), &poly); 

        poly_vec.push_back(poly); 
    }

    //now that we have created the polys, we can create the NPolyModel object
    NPolyModel pmod(poly_vec); 

}