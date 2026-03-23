#include "include/NPolyArray_fit.h"
#include "include/Get_TParameter_from_TFile.h"
#include "include/Add_TParameter_to_TFile.h" 
#include <ROOT/RDataFrame.hxx>
#include <TVector3.h> 
#include <iostream> 

int make_z_hcs_polynomial(int order, const char* path_infile, const char* path_outfile)
{

    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df("tracks_fp", path_infile); 

    std::cout << "Making snapshot..." << std::flush; 
   
    df 
        .Define("z_hcs", [](const TVector3& vtx){ return vtx.z(); }, {"position_vtx"})
        .Snapshot("tracks_fp", "temp.root", {
            "x_fp",
            "y_fp", 
            "dxdz_fp",
            "dydz_fp",

            "z_hcs"
        }); 
    
    //add TParameter to the file 

    bool is_RHRS = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS").value(); 
    
    auto file = new TFile("temp.root", "UPDATE"); 
    Add_TParameter_to_TFile<bool>("is_RHRS", is_RHRS);
    file->Close(); 
    delete file;  

    std::cout << "done." << std::endl; 

    //now, train the polynomial 
    NPolyArray_fit(order, 
        "temp.root",
        path_outfile,  
        {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"}, {"z_hcs"}
    ); 

    return 0; 
}