#include "TROOT.h"
#include <ROOT/RDataFrame.hxx>
#include <TVector3.h> 
#include <TParameter.h>
#include <ApexOptics.h> 
#include "include/Add_TParameter_to_TFile.h" 
#include "include/Get_TParameter_from_TFile.h"

using namespace std; 
using namespace ROOT::VecOps; 

//this covnversion is for monte-carlo data
int convert_mc_to_fpcoords(const char* path_infile, const char* path_outfile, const char* tree_name="tracks_fp") 
{

    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df(tree_name, path_infile); 

    double hrs_momentum = 1104.; 

    //try to get 'is_RHRS' parameter from our input file
    auto param = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS");
    if (!param.has_value()) {
        Error("convert_mc_to_fpcoords", "Unable to extract 'is_RHRS' parameter from TFile: '%s'", path_infile); 
        return -1; 
    }
    const bool is_RHRS = param.value(); 


    cout << "Creating output file..." << flush; 

    auto df_output = df
        .Define("x_sv",     [](TVector3 v){ return v.x(); }, {"position_sieve"})
        .Define("y_sv",     [](TVector3 v){ return v.y(); }, {"position_sieve"})
        .Define("dxdz_sv",  [](TVector3 v){ return v.x()/v.z(); }, {"momentum_sieve"})
        .Define("dydz_sv",  [](TVector3 v){ return v.y()/v.z(); }, {"momentum_sieve"})
        .Define("dpp_sv",   [hrs_momentum](TVector3 v){ return (v.Mag()-hrs_momentum)/hrs_momentum; }, {"momentum_sieve"}) 

        .Define("x_q1",     [](TVector3 v){ return v.x(); }, {"position_Q1"})
        .Define("y_q1",     [](TVector3 v){ return v.y(); }, {"position_Q1"})
        .Define("dxdz_q1",  [](TVector3 v){ return v.x()/v.z(); }, {"momentum_Q1"})
        .Define("dydz_q1",  [](TVector3 v){ return v.y()/v.z(); }, {"momentum_Q1"})
        .Define("dpp_q1",   [hrs_momentum](TVector3 v){ return (v.Mag()-hrs_momentum)/hrs_momentum; }, {"momentum_Q1"})
        
        .Define("position_vtx_scs", [is_RHRS](TVector3 vtx){ return ApexOptics::HCS_to_SCS(is_RHRS, vtx); }, {"position_vtx"})

        .Snapshot(tree_name, path_outfile, {
                "x_fp",
                "y_fp",
                "dxdz_fp",
                "dydz_fp",
                
                "x_q1",
                "y_q1",
                "dxdz_q1",
                "dydz_q1",
                "dpp_q1",

                "x_sv",
                "y_sv",
                "dxdz_sv",
                "dydz_sv",
                "dpp_sv",

                "position_vtx",
                "position_vtx_scs",
                "momentum_vtx"
            });

    cout << "done." << endl; 

    cout << "Adding 'is_RHRS' parameter..." << flush; 

    auto file = new TFile(path_outfile, "UPDATE"); 
    Add_TParameter_to_TFile<bool>("is_RHRS", is_RHRS); 
    file->Close(); 

    cout << "done." << endl; 

    return 0; 
}