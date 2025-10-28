#include <PolynomialCut.h> 
#include <ROOT/RDataFrame.hxx>
#include <TVector3.h>
#include <ApexOptics.h>
#include <ROOT/RVec.hxx>
#include <TStyle.h>
#include <TCanvas.h> 
#include "include/Add_TParameter_to_TFile.h"
#include "include/RDFNodeAccumulator.h"
#include <PolynomialCut.h> 


using ApexOptics::Trajectory_t; 

using namespace std; 

//this is basic script to launch the 'polynomial cut' app, based on cuts in the z_hcs & y_sv space.  
int create_polycut_from_zycut(const char* path_infile) 
{
    const char* const here = "create_polycut_from_zycut"; 

    ROOT::RDataFrame df("tracks_fp", path_infile); 
    RDFNodeAccumulator rna(df); 

    //fetch which arm we're using, and catch any errors that may pop up 
    auto is_RHRS_optional = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS"); 
    
    if (!is_RHRS_optional.has_value()) {
        Error(here, "Unable to fetch 'is_RHRS' param from file."); 
        return -1; 
    }
    const bool is_RHRS = is_RHRS_optional.value(); 

    rna.Define("Xsv", [](double x, double y, double dxdz, double dydz)
        {
            return ApexOptics::RVec_to_Trajectory_t({x, y, dxdz, dydz}); 
        }, {"x_sv", "y_sv", "dxdz_sv", "dydz_sv"});

    rna.Define("Xhcs", [is_RHRS](Trajectory_t Xsv)
        {
            return ApexOptics::SCS_to_HCS(is_RHRS, Xsv); 
        }, {"Xsv"});

    //the new optics runs don't have this variable, so we need to define it for them. 
    rna.DefineIfMissing("y_hcs_corrected", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}); 

    rna.Define("z_reco_horizontal", [is_RHRS](double y_hcs_corrected, Trajectory_t Xhcs)
        {
            return - ( Xhcs.y - y_hcs_corrected ) / Xhcs.dydz; 
        }, {"y_hcs_corrected", "Xhcs"});
    
    rna.Define("z_reco_vertical", [is_RHRS](TVector3 vtx, Trajectory_t Xhcs)
        {
            return - ( Xhcs.x - vtx.x() ) / Xhcs.dxdz; 
        }, {"position_vtx", "Xhcs"}); 


    auto hist_zhcs_vert_dydz  = rna.Get()
        .Histo2D<double>({"h_zv_dydz", "vertical-plane (x_{hcs}=0) z-reco;z_{hcs};dy/dz_{sv}", 200, -0.5,+0.5, 200, -0.04,+0.04}, "z_reco_vertical", "dydz_sv"); 

    auto hist_zhcs_horz_dydz  = rna.Get()
        .Histo2D<double>({"h_zh_dydz", "horizontal-plane (y_{hcs}=0) z-reco;z_{hcs};dx/dz_{sv}", 200, -0.8,+0.7, 200, -0.08,+0.08}, "z_reco_horizontal", "dydz_sv"); 
    

    auto hist_zhcs      = rna.Get()
        .Histo1D<double>({"h_z", "z_{hcs}", 200, -0.5, 0.5}, "z_reco_horizontal"); 
    
    auto hist_x_y       = rna.Get()
        .Histo2D<double>({"h_xy", "x_{sv} vs y_{sv}", 200, -0.05,0.05, 200, -0.05,0.05}, "x_sv", "y_sv");
        
    auto hist_dx_dy      = rna.Get()
        .Histo2D<double>({"h_dxdy", "dx/dz_{sv} vs dy/dz_{sv}", 200, -0.05,0.05, 200, -0.05,0.05}, "dxdz_sv", "dydz_sv");
    
    gStyle->SetOptStat(0);
    
    TCanvas* c=nullptr; 

    PolynomialCut::InteractiveApp((TH2*)hist_zhcs_vert_dydz->Clone("hclone"), "col", kSunset); 
    return 0; 
    
    c = new TCanvas("c1", "z_hcs reconstructions", 1500, 700); 
    c->Divide(2,1); 
    c->cd(1); hist_zhcs_vert_dydz->DrawCopy("col"); 
    c->cd(2); hist_zhcs_horz_dydz->DrawCopy("col"); 

    c = new TCanvas("c2", "xy and angles", 1500, 700); 
    c->Divide(2,1); 
    c->cd(1); hist_x_y->DrawCopy("col"); 
    c->cd(2); hist_dx_dy->DrawCopy("col"); 


    return 0; 
}