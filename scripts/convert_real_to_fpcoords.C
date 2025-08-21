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

    //a very rough cut on sum of all cerenkov channels. 
    const double min_cerenkov_sum = 1e3; 
    string branch_cerenkov_adc = (is_RHRS ? "R" : "L") + "_cer_a_c"; 

    vector<ROOT::RDF::RNode> nodes{ df.Filter([](int n){ return n>0; }, {br_n_tracks}) }; 

    for (auto& branch : fp_coords) { 
        
        const char* br_new = branch.first.c_str(); 

        branch.second = (is_RHRS ? "R" : "L") + branch.second;
        const char* br_old = branch.second.c_str(); 

        nodes.push_back( nodes.back().Define(br_new,   [](const RVec<double>& v){ return v.at(0); }, {br_old}) );

    }

    nodes.back()

        //perform a very basic cut on sum of cerenkov ADCs 
        .Define("cer_sum", [](const ROOT::RVec<double>& cer_a_c)
        {
            double val=0.; 
            for (double adc : cer_a_c) if (fabs(adc) < 5e4) val += adc; 
            return val; 
        }, {branch_cerenkov_adc.c_str()})

        .Filter([min_cerenkov_sum](double cer_sum){ return cer_cum > min_cerenkov_sum; }, {"cer_sum"})

        //perform the conversion from 'transport' coordinates (tra) to 'focal-plane' coordinates (fp)
        .Redefine("dxdz_fp", [](double x_tra, double dxdz_tra){
            return dxdz_tra - x_tra/6.; 
        }, {"x_fp", "dxdz_fp"})

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

    auto param_is_RHRS     = new TParameter<bool>  ("is_RHRS", is_RHRS);      
    auto param_cer_sum_cut = new TParameter<double>("min_cerenkov_sum", min_cerenkov_sum); 
        
    param_is_RHRS    ->Write();
    param_cer_sum_cut->Write(); 

    file->Close(); 

    cout << "done." << endl; 

    return 0; 
}