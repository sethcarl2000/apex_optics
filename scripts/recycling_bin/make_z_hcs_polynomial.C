#include <ApexOptics.h>
#include "include/NPolyArray_fit.h"
#include "include/Get_TParameter_from_TFile.h"
#include "include/Add_TParameter_to_TFile.h" 
#include <ROOT/RDataFrame.hxx>
#include <TVector3.h> 
#include <iostream> 

int make_z_hcs_polynomial(int order, const char* path_infile, const char* path_outfile)
{
    using ApexOptics::Trajectory_t; 

    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df("tracks_fp", path_infile); 

    const bool is_RHRS = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS").value(); 

    std::cout << "Making snapshot..." << std::flush; 
   
    df 
        .Define("z_hcs_reco", [is_RHRS](double dxdz, double dydz, double x, double y, TVector3 vtx_hcs)
        { 
            Trajectory_t Xsv{ x, y, dxdz, dydz }; 

            auto Xhcs = ApexOptics::SCS_to_HCS(is_RHRS, Xsv); 

            return - ( Xhcs.x ) / Xhcs.dxdz; 

        }, {"dxdz_sv", "dydz_sv", "x_sv", "y_sv", "position_vtx"})

        .Define("z_hcs", [](const TVector3& vtx){ return vtx.z(); }, {"position_vtx"})

        .Snapshot("tracks_fp", "temp.root", {
            "x_fp",
            "y_fp", 
            "dxdz_fp",
            "dydz_fp",
            "dpp_sv",

            "z_hcs_reco",
            "z_hcs"
        }); 
    
    auto file = new TFile("temp.root", "UPDATE"); 
    Add_TParameter_to_TFile<bool>("is_RHRS", is_RHRS);
    file->Close(); 
    delete file;  

    std::cout << "done." << std::endl; 

    //now, train the polynomial 
    NPolyArray_fit(order, 
        "temp.root",
        path_outfile,  
        {"x_fp", "y_fp", "dxdz_fp", "dydz_fp", "dpp_sv"}, {"z_hcs_reco"}
    ); 

    return 0; 
}