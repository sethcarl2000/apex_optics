#include "TROOT.h"
#include <ROOT/RDataFrame.hxx>
#include <string> 

using namespace std; 
using namespace ROOT::VecOps; 

//this covnversion is for real optics data data
int convert_real_to_fpcoords(const bool is_RHRS, const char* path_infile, const char* path_outfile) 
{
    const char* tree_name = "track_data"; 

    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df(tree_name, path_infile); 

    cout << "Creating output file..." << flush; 

    vector<pair<string, string>> fp_coords = {
        {"x_fp",    "_tr_tra_x" },
        {"y_fp",    "_tr_tra_y" },
        {"dxdz_fp", "_tr_tra_th"},
        {"dydz_fp", "_tr_tra_ph"},
    }; 

    const char* br_n_tracks = (is_RHRS ? "R_tr_n" : "L_tr_n"); 


    vector<ROOT::RDF::RNode> nodes{ df.Filter([](int n){ return n>0; }, {br_n_tracks}) }; 

    for (auto& branch : fp_coords) { 
        
        const char* br_new = branch.first.c_str(); 

        branch.second = (is_RHRS ? "R" : "L") + branch.second;
        const char* br_old = branch.second.c_str(); 

        nodes.push_back( nodes.back().Define(br_new,   [](const RVec<double>& v){ return v.at(0); }, {br_old}) );

    }

    nodes.back()
        .Define("position_vtx", [](TVector3 vtx){ return vtx; }, {"react_vertex"})
        .Snapshot("tracks_fp", path_outfile, {
            "x_fp",
            "y_fp",
            "dxdz_fp",
            "dydz_fp",

            "position_vtx"
        });

    cout << "done." << endl; 

    cout << "Adding 'is_RHRS' parameter..." << flush; 

  
    auto file = new TFile(path_outfile, "UPDATE"); 

    auto param_is_RHRS = new TParameter<bool>("is_RHRS", is_RHRS); 
    param_is_RHRS->Write();
    file->Close(); 

    cout << "done." << endl; 

    return 0; 
}