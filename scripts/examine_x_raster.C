#include <PolynomialCut.h> 
#include <ROOT/RDataFrame.hxx>
#include <TVector3.h>
#include <TPad.h> 
#include <ApexOptics.h>
#include <ROOT/RVec.hxx>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>  
#include "include/Add_TParameter_to_TFile.h"
#include "include/RDFNodeAccumulator.h"
#include <PolynomialCut.h> 


using ApexOptics::Trajectory_t; 

using namespace std; 

//this is basic script to launch the 'polynomial cut' app, based on cuts in the z_hcs & y_sv space.  
int examine_x_raster(const char* path_infile) 
{
    const char* const here = "examine_x_raster"; 

    ROOT::RDataFrame df("tracks_fp", path_infile); 
    RDFNodeAccumulator rna(df); 

    //fetch which arm we're using, and catch any errors that may pop up 
    auto is_RHRS_optional = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS"); 
    
    if (!is_RHRS_optional.has_value()) {
        Error(here, "Unable to fetch 'is_RHRS' param from file."); 
        return -1; 
    }
    const bool is_RHRS = is_RHRS_optional.value(); 

    PolynomialCut pcut_V2, pcut_H3, pcut_V3; 
    try { 
        cout << "Parsing polynomial cuts..." << flush; 

        pcut_V2.Parse_dbfile("data/polycuts/V2-4775.dat");
        pcut_H3.Parse_dbfile("data/polycuts/H3-4775.dat");
        pcut_V3.Parse_dbfile("data/polycuts/V3-4775.dat");
    } 
    catch (const PolynomialCut::DBFileException& e) {
        Error(here, "Problem parsing cutfile. Exception:\n %s", e.what());
        return -1;  
    }
    cout << "done." << endl; 



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

    
    auto df_V2 = rna.Get().Filter([&pcut_V2](double z_hcs, double dydz){ return pcut_V2.IsInside(z_hcs,dydz); }, {"z_reco_vertical", "dydz_sv"}); 
    auto df_H3 = rna.Get().Filter([&pcut_H3](double z_hcs, double dydz){ return pcut_H3.IsInside(z_hcs,dydz); }, {"z_reco_vertical", "dydz_sv"}); 
    auto df_V3 = rna.Get().Filter([&pcut_V3](double z_hcs, double dydz){ return pcut_V3.IsInside(z_hcs,dydz); }, {"z_reco_vertical", "dydz_sv"}); 
    
    
    const double range[] = { 27.5e3, 32.0e3 }; 
    const int nbins = 100; 


    auto hist_V2 = df_V2
        .Histo1D<double>({"h_V2", "V2;Raster2 current (x);", nbins, range[0], range[1]}, "Raster2_current_x");
        
    auto hist_H3 = df_H3
        .Histo1D<double>({"h_H3", "H3;Raster2 current (x);", nbins, range[0], range[1]}, "Raster2_current_x");
    
    auto hist_V3 = df_V3
        .Histo1D<double>({"h_V3", "V3;Raster2 current (x);", nbins, range[0], range[1]}, "Raster2_current_x");
        
    auto hist_all = rna.Get()
        .Histo1D<double>({"h_all", "all;Raster2 current (x);", nbins, range[0], range[1]}, "Raster2_current_x");
        

    gStyle->SetOptStat(0);
    
    TCanvas* c=nullptr; 

    c = new TCanvas("c1", "z_hcs reconstructions", 700, 700); 
    
    hist_all->SetLineStyle(kDashed); 
    hist_all->SetFillStyle(0);
    hist_all->SetLineColor(kBlack); 
    hist_all->GetXaxis()->SetNdivisions(-5);
    hist_all->SetTitle("All;Raster-2 raw X-current ADC (arb. units);"); 
    //hist_all->DrawCopy(); 
    
    hist_V2->SetLineColor(kRed); 
    hist_V2->SetFillStyle(3003); 
    hist_V2->SetFillColor(kRed);
    hist_V2->SetTitle("V2;Raster-2 raw X-current ADC (arb. units);"); 
    hist_V2->DrawCopy(); 

    hist_H3->SetLineColor(kBlack); 
    hist_H3->SetLineStyle(kDotted);
    hist_H3->SetFillStyle(3005); 
    hist_H3->SetFillColor(kBlack);  
    
    //hist_H3->DrawCopy("SAME"); 
    
    hist_V3->SetLineColor(kBlue);
    hist_V3->SetFillStyle(3004); 
    hist_V3->SetFillColor(kBlue); 
    hist_V3->DrawCopy("SAME"); 

    
    gPad->BuildLegend(); 

    return 0; 
}