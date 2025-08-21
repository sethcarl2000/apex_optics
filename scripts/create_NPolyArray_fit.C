#include <memory>
#include <string> 
#include <map>
#include <fstream>
#include <iomanip> 
#include "TFile.h"
#include "TVector3.h"
#include "TParameter.h"


using namespace std; 
   
struct Track_t {
    double x,y,dxdz,dydz,dpp; 
};

//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
//if you want to avoid creating (or overwriting) a db-file, just enter "" as the db-file name, and only the histograms will be drawn instead. 
int create_NPolyArray_fit(  const int poly_order=2,
                            const char* path_infile="",
                            const char* stem_outfile="",  
                            const vector<string> inputs = {},
                            const vector<string> outputs = {}, 
                            const char* tree_name="tracks_fp" ) 
{
    const char* const here = "fit_points_mc_forward"; 

    if (inputs.empty() || outputs.empty()) {
        Error(here, "inputs and/or outputs are empty");
        return -1; 
    }
    
    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(here, "libApexOptics could not be loaded."); 
        return 1; 
    }

    auto infile = new TFile(path_infile, "READ");

    //check if the infile could be opened
    if (!infile || infile->IsZombie()) {
        Error(here, "root file '%s' could not be opened.", path_infile); 
        return 1; 
    }

    //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
    TParameter<bool>* param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
    if (!param_is_RHRS) {
        Error(here, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
        return 1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 

    infile->Close(); 
    delete infile; 

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame

    auto df_ptr = unique_ptr<ROOT::RDataFrame>(nullptr); 
    try { df_ptr = unique_ptr<ROOT::RDataFrame>(new ROOT::RDataFrame(tree_name, path_infile)); }
    catch (const std::invalid_argument& e) {
        Error(here, "Trying to create RDataFrame threw std::invalid_argument exception.\n what(): %s", e.what()); 
        return -1; 
    } 
    ROOT::RDataFrame& df = *df_ptr; 
    
    const size_t DoF_in  = inputs.size();  
    const size_t DoF_out = outputs.size();  

    cout << "Creating polynomials for inputs => outputs..." << flush; 
    
    //for each of the branches in the 'branches_sv' created above, make a polynomial which takes all the branches
    // of the 'branches_fp' vec 
    map<string,NPoly*> polymap = ApexOptics::Create_NPoly_fit(df, poly_order, inputs, outputs);
    
    cout << "done." << endl;    

    string path_outfile(stem_outfile);

    //create the output file ___________________________________________
    if (path_outfile != "") {
        cout << "Creating dbfiles for polynomials..." << flush; 
        
        ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap); 

        cout << "done." << endl; 

    } else { 

        cout << "Skipping db-file creation." << endl; 
    }
    //____________________________________________________________________________

    //delete our poly models
    for (auto it = polymap.begin(); it != polymap.end(); it++ ) delete it->second;
    
    return 0;
}