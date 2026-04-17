
#ifndef test_opticsModel_C
#define test_opticsModel_C 

//APEX headers
#include "include/RDFNodeAccumulator.h"
#include "include/Get_TParameter_from_TFile.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include <ApexOptics.h> 
#include <ModularOpticsModel.h>
//ROOT headers
#include <TParameter.h>
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h> 
//stdlib headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>

using ApexOptics::Trajectory_t; 

/// @brief draws graphs for the given optics model 
/// @param model model to be tested
int test_opticsModel(
    ModularOpticsModel* model, 
    std::string path_infile, 
    std::string target=""
) 
{
    if (!model) {
        Error(__func__, "Passed null optics model ptr"); 
        return -1; 
    }

    const char* tree_name = "track_data";

    using namespace std; 
    using namespace ROOT::VecOps; 

    auto infile = new TFile(path_infile.data(), "READ");

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(__func__, "libApexOptics could not be loaded."); 
        return 1; 
    }

    //const bool is_RHRS = Get_TParameter_from_TFile<bool>(path_infile.data(), "is_RHRS").value(); 
    

    ROOT::EnableImplicitMT(); 
    Info(__func__, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //this is a little helper class which is meant to avoid some of the awkward syntax typically associated with RDataFrames creation. 
    auto rna = RDFNodeAccumulator(df);  

    //define output branches
    rna = model->DefineOutputs(rna.Get()); 

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(__func__, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    //create both histograms
    auto R_hist_xy = rna.Get()
        .Histo2D({"h_xy",     "(RHRS) Sieve-plane projection;x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, +0.025}, "R_x_sv", "R_y_sv");
    
    auto R_hist_angles = rna.Get()
        .Histo2D({"h_angles", "(RHRS) Sieve-plane projection;dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, +0.04}, "R_dxdz_sv", "R_dydz_sv"); 
    

    //create both histograms
    auto L_hist_xy = rna.Get()
        .Histo2D({"h_xy",     "(LHRS) Sieve-plane projection;x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, +0.025}, "L_x_sv", "L_y_sv");
    
    auto L_hist_angles = rna.Get()
        .Histo2D({"h_angles", "(LHRS) Sieve-plane projection;dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, +0.04}, "L_dxdz_sv", "L_dydz_sv"); 
    
    auto his_xz_hcs = rna.Get()
        .Define("x_hcs", [](TVector3 vtx){return vtx.x();}, {"reco_position_vtx"})
        .Define("z_hcs", [](TVector3 vtx){return vtx.z();}, {"reco_position_vtx"})
        .Histo2D({"h_xz", "React vertex reconstruction (x & z);z_{hcs} (m);x_{hcs} (m)", 200, -0.30,+0.30, 200,-0.30,+0.30}, "z_hcs", "x_hcs"); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 

    //PolynomialCut::InteractiveApp((TH2*)hist_z_y->Clone("hclone"), "col2", kSunset); 
    //return 0; 

    TCanvas* c; 
    
    c = new TCanvas("c1_R", Form("data: %s",path_infile.data()), 1200, 600); 
    gPad->SetLeftMargin(0.15); 
    c->Divide(2,1, 0.01,0.01); 

    (c->cd(1))->SetLeftMargin(0.15); R_hist_xy->DrawCopy("col2"); 
    (c->cd(2))->SetLeftMargin(0.15); R_hist_angles->DrawCopy("col2");


    c = new TCanvas("c1_L", Form("data: %s",path_infile.data()), 1200, 600); 
    gPad->SetLeftMargin(0.15); 
    c->Divide(2,1, 0.01,0.01); 

    (c->cd(1))->SetLeftMargin(0.15); L_hist_xy->DrawCopy("col2"); 
    (c->cd(2))->SetLeftMargin(0.15); L_hist_angles->DrawCopy("col2");


    c = new TCanvas("c1_xz", Form("data: %s",path_infile.data()), 1200, 600); 
    gPad->SetLeftMargin(0.15); 

    his_xz_hcs->DrawCopy("col"); 

    return 0; 
}

#endif 