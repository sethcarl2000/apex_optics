#include <memory>
#include <string> 
#include <map>
#include <fstream>
#include <iomanip> 
#include "TFile.h"
#include "TVector3.h"
#include "TParameter.h"


using namespace std; 


//creates db '.dat' files for two sepearate polynomials:
// the fp_q1 polynomials map from FOCAL PLANE coordinates to Q1 FRONT coordinates.
// the q1_sv polynomials map from Q1 FRONT coordinates to SIEVE coordinates. 
int fitpoints_mc_fp_q1_sv(  const int poly_fpq1_order=2,
                            const int poly_q1sv_order=2,
                            const char* path_infile="",
                            const char* stem_outfile="data/csv/db_mc",  
                            const char* tree_name="tracks_fp" ) 
{
    const char* const here = "fitpoints_mc_fp_q1_sv"; 

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
    //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
    TParameter<bool>* param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
    if (!param_is_RHRS) {
        Error(here, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
        return 1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 

    delete tree; 
    delete is_RHRS; 
    infile->Close(); 
    delete infile; 

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    
    const int hrs_momentum = 1104.0; 

    auto df_output = df

        .Define("x_q1",      [](TVector3 v){ return v.x(); },        {"position_Q1"})
        .Define("y_q1",      [](TVector3 v){ return v.y(); },        {"position_Q1"})
        .Define("dxdz_q1",   [](TVector3 v){ return v.x()/v.z(); },  {"momentum_Q1"})
        .Define("dydz_q1",   [](TVector3 v){ return v.y()/v.z(); },  {"momentum_Q1"})
        .Define("dpp_q1",    [hrs_momentum](TVector3 v){ return (v.Mag()-hrs_momentum)/hrs_momentum; }, {"momentum_Q1"})

        .Define("x_sv",      [](TVector3 v){ return v.x(); },        {"position_sieve"})
        .Define("y_sv",      [](TVector3 v){ return v.y(); },        {"position_sieve"})
        .Define("dxdz_sv",   [](TVector3 v){ return v.x()/v.z(); },  {"momentum_sieve"})
        .Define("dydz_sv",   [](TVector3 v){ return v.y()/v.z(); },  {"momentum_sieve"}); 
         
    //now, put these outputs into a vector (so we know to make a seperate polynomial for each of them). 
    vector<string> branches_fp = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 
    
    vector<string> branches_q1 = {
        "x_q1", 
        "y_q1",
        "dxdz_q1",
        "dydz_q1",
        "dpp_q1"
    };  

    vector<string> branches_sv = {
        "x_sv", 
        "y_sv",
        "dxdz_sv",
        "dydz_sv" 
    }; 


    //these 'result pointers' will let us see the result for each element of the least-squares fit matrix
    //first, find the best coffeficients for the polynomials.
    
    cout << "Creating polynomials for fp => q1..." << endl; 
    
    map<string,NPoly*> polymap_fpq1;     
    for (const string& output : branches_q1 ) 
        polymap_fpq1[output] = ApexOptics::Create_NPoly_fit(df_output, poly_fpq1_order, branches_fp, output.data());
    
    cout << "done." << endl;    

    
    cout << "Creating polynomials for q1 => sv..." << endl; 

    map<string,NPoly*> polymap_q1sv;     
    for (const string& output : branches_sv ) 
        polymap_q1sv[output] = ApexOptics::Create_NPoly_fit(df_output, poly_q1sv_order, branches_q1, output.data());
    
    cout << "done." << endl;  


    //construct the full path to the outfile from the stem provided
    string path_outfile; 
    char buffer[25];
    
    //create the fp => q1 output file ___________________________________________
    path_outfile = string(stem_outfile); 

    //specify the arm to use 
    path_outfile += "_fp_q1"; 

    path_outfile += (is_RHRS ? "_R" : "_L"); 

    //specify the order of the polynomial
    sprintf(buffer, "_%iord", poly_fpq1_order); 
    path_outfile += string(buffer); 

    path_outfile += ".dat";

    //write all elements of the output file
    ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap_fpq1); 


    //create the q1 => sv output file ___________________________________________
    path_outfile = string(stem_outfile); 

    //specify the arm to use 
    path_outfile += "_q1_sv"; 

    path_outfile += (is_RHRS ? "_R" : "_L"); 

    //specify the order of the polynomial
    sprintf(buffer, "_%iord", poly_q1sv_order); 
    path_outfile += string(buffer); 

    path_outfile += ".dat";

    //write all elements of the output file
    ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap_q1sv); 

    //delete our template polynomial, and all model polynomials
    for (auto it = polymap_fpq1.begin(); it != polymap_fpq1.end(); it++ ) delete it->second;

    for (auto it = polymap_q1sv.begin(); it != polymap_q1sv.end(); it++ ) delete it->second;


    return 0;
}