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
int fitpoints_mc_fp_sv( const int poly_order=2,
                        const char* path_infile="",
                        const char* stem_outfile="",  
                        const char* tree_name="tracks_fp" ) 
{
    const char* const here = "fit_points_mc_forward"; 

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

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    
    const int nDoF = 4; //the 4 DoF are: x_fp, y_fp, dxdz_fp, dydz_fp

    auto poly = new NPoly(nDoF,poly_order); 

    //get number of coefficients
    const int n_elems = poly->Get_nElems(); 

    const double hrs_momentum = 1104.0; 

    //this will be used for the least-squares calculation to find the best coefficients
    RMatrix A_init(n_elems, n_elems, 0.); 

    //Define all of the branches we want to create models to map between
    auto df_output = df

        .Define("x_sv",      [](TVector3 v){ return v.x(); },        {"position_sieve"})
        .Define("y_sv",      [](TVector3 v){ return v.y(); },        {"position_sieve"})
        .Define("dxdz_sv",   [](TVector3 v){ return v.x()/v.z(); },  {"momentum_sieve"})
        .Define("dydz_sv",   [](TVector3 v){ return v.y()/v.z(); },  {"momentum_sieve"})
        .Define("dpp_sv",    [hrs_momentum](TVector3 v){ return (v.z()-hrs_momentum)/hrs_momentum; }, {"momentum_sieve"});
    
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

    cout << "Creating polynomials for fp => sv..." << flush; 
    
    //for each of the branches in the 'branches_sv' created above, make a polynomial which takes all the branches
    // of the 'branches_fp' vec 
    map<string,NPoly*> polymap =       
        ApexOptics::Create_NPoly_fit(df_output, poly_order, branches_fp, branches_sv);
    
    cout << "done." << endl;    

    string path_outfile = "";

    //create the fp => sv output file ___________________________________________
    if (string(stem_outfile) != "") {
        cout << "Creating dbfiles for polynomials..." << flush; 
        
        path_outfile = string(stem_outfile); 
        char buffer[50]; 

        //specify the arm to use 
        path_outfile += "_fp_sv"; 

        path_outfile += (is_RHRS ? "_R" : "_L"); 

        //specify the order of the polynomial
        sprintf(buffer, "_%iord", poly_order); 
        path_outfile += string(buffer); 

        path_outfile += ".dat";
    
        ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap); 

        cout << "done." << endl; 

    } else { cout << "Skipping db-file creation." << endl; }
    //____________________________________________________________________________

    
    //draw the reults of all models
    //reconstruct the sieve coordinates
    NPolyArray* parr = new NPolyArray({
        polymap["x_sv"],
        polymap["y_sv"],
        polymap["dxdz_sv"],
        polymap["dydz_sv"],
        polymap["dpp_sv"]
    }); 
    
    auto df_sieve_reco = df_output 
        
        //Define the reconstructed sieve coordinates using our model
        .Define("reco_X_sv",    [parr](double x, double y, double dxdz, double dydz)
        {
            auto X_sv = parr->Eval({x, y, dxdz, dydz}); 
            return Track_t{ 
                .x      =X_sv[0], 
                .y      =X_sv[1], 
                .dxdz   =X_sv[2],
                .dydz   =X_sv[3],
                .dpp    =X_sv[4]
            };
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"}); 
    
        
    vector<ROOT::RDF::RNode> error_nodes{ df_sieve_reco }; 
    
    auto Define_Track_t_branch = [&error_nodes](const char* coord_name, double Track_t::*coord)
    {   
        auto new_node = error_nodes.back() 
            .Define(coord_name, [coord](const Track_t& trk) { return trk.*coord; }, {"reco_X_sv"});

        error_nodes.push_back(new_node); 
        return; 
    };
    
    Define_Track_t_branch("reco_x_sv",      &Track_t::x);
    Define_Track_t_branch("reco_y_sv",      &Track_t::y);
    Define_Track_t_branch("reco_dxdz_sv",   &Track_t::dxdz);
    Define_Track_t_branch("reco_dydz_sv",   &Track_t::dydz);
    Define_Track_t_branch("reco_dpp_sv",    &Track_t::dpp);
    

    for (auto it = polymap.begin(); it != polymap.end(); it++) {

        //get the polynomial and its name
        const char* poly_name = it->first.data(); 
        NPoly *poly           = it->second; 
        
        //name the new output branch 
        char br_output_name[50];
        sprintf(br_output_name, "error_%s", poly_name);

        auto new_node = error_nodes.back() 

            .Define(br_output_name, [poly](double target, double x, double y, double dxdz, double dydz)
            {
                return (target - poly->Eval({x, y, dxdz, dydz}))*1e3; 

            }, {poly_name, "x_fp", "y_fp", "dxdz_fp", "dydz_fp"}); 

        error_nodes.push_back(new_node); 
    }

    auto df_error = error_nodes.back(); 

    auto h_x    = df_error.Histo1D({"h_x", "Error of x_sv;mm", 200, -10, 10}, "error_x_sv"); 
    auto h_y    = df_error.Histo1D({"h_y", "Error of y_sv;mm", 200, -10, 10}, "error_y_sv"); 
    auto h_dxdz = df_error.Histo1D({"h_dxdz", "Error of dxdz_sv;mrad", 200, -2, 2}, "error_dxdz_sv"); 
    auto h_dydz = df_error.Histo1D({"h_dydz", "Error of dydz_sv;mrad", 200, -2, 2}, "error_dydz_sv"); 
    
    auto h_xy_sieve = df_error.Histo2D<double>({"h_xy_sieve", "Reconstructed sieve-coords;x_sv;y_sv", 250, -45e-3, 45e-3, 250, -45e-3, 45e-3}, "reco_x_sv", "reco_y_sv"); 
    

    char b_c_title[120]; 
    sprintf(b_c_title, "Errors of different coords: %s", path_outfile.data()); 
    auto c = new TCanvas("c", b_c_title, 1200, 800); 

    c->Divide(2,2, 0.005,0.005); 
    
    c->cd(1); 
    h_x->DrawCopy(); 
    c->cd(2); 
    h_y->DrawCopy(); 
    c->cd(3); 
    h_dxdz->DrawCopy(); 
    c->cd(4); 
    h_dydz->DrawCopy(); 


    new TCanvas("c2", b_c_title); 

    gStyle->SetPalette(kSunset); 
    h_xy_sieve->SetStats(0); 

    h_xy_sieve->DrawCopy("col2"); 

    //delete our template polynomial
    delete poly;
    delete parr;  

    //delete our poly models
    for (auto it = polymap.begin(); it != polymap.end(); it++ ) delete it->second;

    return 0;
}