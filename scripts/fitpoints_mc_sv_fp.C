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
int fitpoints_mc_sv_fp( bool is_RHRS=false,
                        const int poly_order=2,
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
/*
        .Define("X_elems", [poly](double x, double y, double dxdz, double dydz)
    {
        return poly->Eval_noCoeff({x, y, dxdz - x/6., dydz}); 
    }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"})
*/
        //this is the only difference between VDC TRANSPORT COORDINATES (tra) and FOCAL PLANE COORDINATES (fp)
        .Redefine("dxdz_fp", [](double x_tra, double dxdz_tra){return dxdz_tra - x_tra/6.;}, {"x_fp", "dxdz_fp"})

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

    cout << "Creating polynomials for fp => sv..." << endl; 
    
    //for each of the branches in the 'branches_sv' created above, make a polynomial which takes all the branches
    // of the 'branches_fp' vec 
    map<string,NPoly*> polymap;     
    for (const string& output : branches_fp ) {
        polymap[output] = ApexOptics::Create_NPoly_fit( df_output, //the dataframe object with all of our branches defined 
                                                        poly_order, //the order of the NPoly to create 
                                                        branches_sv, //the input branches 
                                                        output.data()); //the output branch for this NPoly to target
    
                                                    }

    cout << "done." << endl;    


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
    vector<ROOT::RDF::RNode> error_nodes{ df_output }; 
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
                return (target - poly->Eval({x, y, dxdz, dydz, dpp}))*1e3; 

            }, {poly_name, "x_sv", "y_sv", "dxdz_sv", "dydz_sv", "dpp_sv"}); 

        error_nodes.push_back(new_node); 
    }

    auto df_error = error_nodes.at(error_nodes.size()-1); 

    char b_c_title[120]; 
    sprintf(b_c_title, "Errors of different coords: %s", path_outfile.data()); 
    auto c = new TCanvas("c", b_c_title, 1200, 800); 

    c->Divide(2,2, 0.005,0.005); 
    
    c->cd(1); 
    auto h_x = df_error.Histo1D({"h_x", "Error of x_fp;mm", 200, -10, 10}, "error_x_fp"); 
    h_x->DrawCopy(); 
    c->cd(2); 
    auto h_y = df_error.Histo1D({"h_y", "Error of y_fp;mm", 200, -10, 10}, "error_y_fp"); 
    h_y->DrawCopy(); 
    c->cd(3); 
    auto h_dxdz = df_error.Histo1D({"h_dxdz", "Error of dxdz_fp;mrad", 200, -2, 2}, "error_dxdz_fp"); 
    h_dxdz->DrawCopy(); 
    c->cd(4); 
    auto h_dydz = df_error.Histo1D({"h_dydz", "Error of dydz_fp;mrad", 200, -2, 2}, "error_dydz_fp"); 
    h_dydz->DrawCopy(); 


    //delete our template polynomial
    delete poly; 

    //delete our poly models
    for (auto it = polymap.begin(); it != polymap.end(); it++ ) delete it->second;

    return 0;
}