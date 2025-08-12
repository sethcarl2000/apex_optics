#include <memory>
#include <string> 
#include <map>
#include <fstream>
#include <iomanip> 
#include "TFile.h"
#include "TVector3.h"
#include "TParameter.h"


using namespace std; 

//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
int fitpoints_mc_sv_fp(  const int poly_order=2,
                         const char* path_infile="",
                         const char* stem_outfile="data/csv/db_mc",  
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

    //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
    TParameter<bool>* param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
    if (!param_is_RHRS) {
        Error(here, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
        return 1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 

    delete tree; 
    infile->Close(); 
    delete infile; 

    const double hrs_momentum = 1104.0; 

    //now, put these outputs into a vector (so we know to make a seperate polynomial for each of them). 
    vector<string> branches_sv = {
        "x_sv", 
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    };  

    vector<string> branches_fp = {
        "x_fp", 
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    };

    cout << "Creating polynomials for fp => sv..." << endl; 
    
    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //for each of the branches in the 'branches_sv' created above, make a polynomial which takes all the branches
    // of the 'branches_fp' vec 
    map<string,NPoly*> polymap = 
        ApexOptics::Create_NPoly_fit( df, //the dataframe object with all of our branches defined 
                                      poly_order, //the order of the NPoly to create 
                                      branches_sv, //the input branches 
                                      branches_fp ); //the output branch for this NPoly to target

    cout << "done." << endl;    

    NPoly *x_fp = polymap["x_fp"]; 
    

    //create the fp => sv output file ___________________________________________
    string path_outfile = string(stem_outfile); 
    char buffer[50]; 

    //specify the arm to use 
    path_outfile += "_sv_fp"; 

    path_outfile += (is_RHRS ? "_R" : "_L"); 

    //specify the order of the polynomial
    sprintf(buffer, "_%iord", poly_order); 
    path_outfile += string(buffer); 

    path_outfile += ".dat";

    ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap);     
    
    //draw the reults of all models
    vector<ROOT::RDF::RNode> error_nodes{ df }; 
    for (auto it = polymap.begin(); it != polymap.end(); it++) {

        //get the polynomial and its name
        const char* poly_name = it->first.data(); 
        NPoly *poly           = it->second; 
        
        //name the new output branch 
        char br_output_name[50];
        sprintf(br_output_name, "error_%s", poly_name);

        //evaluate the difference between the coordinate (target) and the polynomial meant to model it (poly)
        auto new_node = error_nodes.at(error_nodes.size()-1) 

            .Define(br_output_name, [poly](double target, double x, double y, double dxdz, double dydz, double dpp)
            {
                return ( target - poly->Eval({x, y, dxdz, dydz, dpp}) )*1e3; 

            }, {poly_name, "x_sv", "y_sv", "dxdz_sv", "dydz_sv", "dpp_sv"}); 

        error_nodes.push_back(new_node); 
    }

    auto df_error = error_nodes.at(error_nodes.size()-1); 

    auto h_x    = df_error.Histo1D({"h_x",    "Error of x_fp;mm", 200, -5, 5}, "error_x_fp"); 
    auto h_y    = df_error.Histo1D({"h_y",    "Error of y_fp;mm", 200, -5, 5}, "error_y_fp"); 
    auto h_dxdz = df_error.Histo1D({"h_dxdz", "Error of dxdz_fp;mrad", 200, -5, 5}, "error_dxdz_fp"); 
    auto h_dydz = df_error.Histo1D({"h_dydz", "Error of dydz_fp;mrad", 200, -5, 5}, "error_dydz_fp"); 
    

    char b_c_title[120]; 
    sprintf(b_c_title, "Errors of different coords: %s", path_outfile.data()); 
    auto c = new TCanvas("c", b_c_title, 1200, 800); 

    c->Divide(2,2, 0.005,0.005); 
    
    c->cd(1); h_x->DrawCopy(); 
    c->cd(2); h_y->DrawCopy(); 
    c->cd(3); h_dxdz->DrawCopy(); 
    c->cd(4); h_dydz->DrawCopy(); 


    //delete our poly models
    for (auto it = polymap.begin(); it != polymap.end(); it++ ) delete it->second;

    return 0;
}