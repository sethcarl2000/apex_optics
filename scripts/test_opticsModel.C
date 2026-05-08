
#ifndef test_opticsModel_C
#define test_opticsModel_C 

//APEX headers
#include "include/RDFNodeAccumulator.h"
#include "include/Get_TParameter_from_TFile.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include <ApexOptics.h> 
#include <ModularOpticsModel.h>
//ROOT headers
#include <ROOT/RVec.hxx>
#include <TParameter.h>
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h> 
#include <TGraph.h> 
#include <ROOT/RDataFrame.hxx>
//stdlib headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>  
#include <map>

using ApexOptics::Trajectory_t; 

namespace {
    
    //grab a vector of values from the meta_data tree
    template<typename T> std::vector<T> grab_vector_from_node(ROOT::RNode node, std::string branch_name) {
        return (std::vector<T>)*node.Take<T>(branch_name); 
    }; 
}

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
    using ApexOptics::SCS_to_HCS; 

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

    constexpr double hrs_momentum = 1104.;
    
    //define invariant mass
    rna.Define("reco_invariant_mass", [hrs_momentum](const Trajectory_t& R_Xsv_scs, const Trajectory_t& L_Xsv_scs)
    {
        auto R_Xsv = SCS_to_HCS(true,  R_Xsv_scs);
        auto L_Xsv = SCS_to_HCS(false, L_Xsv_scs);

        TVector3 R_p( R_Xsv.dxdz, R_Xsv.dydz, 1. ); R_p = R_p.Unit() * hrs_momentum*( 1. + R_Xsv.dpp );
        TVector3 L_p( L_Xsv.dxdz, L_Xsv.dydz, 1. ); L_p = L_p.Unit() * hrs_momentum*( 1. + L_Xsv.dpp );
        
        //we're neglecting the electron mass
        double m2 = 2.*( R_p.Mag()*L_p.Mag() - (R_p * L_p) );

        return std::sqrt(m2);

    }, {"R_Xsv_reco", "L_Xsv_reco"}); 

    rna.Define("x_hcs", [](TVector3 vtx){return vtx.x();}, {"reco_position_vtx"});
    rna.Define("z_hcs", [](TVector3 vtx){return vtx.z();}, {"reco_position_vtx"});

    rna = rna.Get().Filter([](double x_hcs, double z_hcs){
        if (z_hcs > 0.300 || z_hcs < -0.300) return false; 
        if (std::fabs(x_hcs) > 0.080) return false; 
        return true; 
    }, {"x_hcs", "z_hcs"}); 
     

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
        .Histo2D({"h_xz", "React vertex reconstruction (x & z);z_{hcs} (m);x_{hcs} (m)", 200, -0.30,+0.30, 200,-0.30,+0.30}, "z_hcs", "x_hcs"); 

    auto hist_m = rna.Get() 
        .Histo1D<double>({"h_m", "Invariant mass (MeV);m_{A} (MeV);counts / MeV", 200, 100, 300}, "reco_invariant_mass"); 


    auto hist_dt = rna.Get()
        .Histo1D({"h_dt", "(T[S2,R] - T[S2,L]) / #sigma_{R-L};", 200, -1, -1}, "L_track_dt"); 
    //set the color pallete
    gStyle->SetPalette(kSunset); 
    //gStyle->SetOptStat(0); 

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

    //look at some meta-data 
    ROOT::RDataFrame df("meta_data", path_infile); 


    /*vector<int> v_run = grab_vector_from_node<int>(df, "run_number"); 
    auto v_n_total = grab_vector_from_node<int>(df, "n_events_total"); 
    auto v_n_coinc = grab_vector_from_node<ULong64_t>(df, "n_events_coinc"); 
    auto v_n_events_track_R = grab_vector_from_node<ULong64_t>(df, "n_events_1raw_R"); */


    c = new TCanvas("c1_xz", Form("data: %s",path_infile.data()), 1200, 600); 
    gPad->SetLeftMargin(0.15); 

    his_xz_hcs->DrawCopy("col"); 

    c = new TCanvas("c_m", Form("data: %s",path_infile.data()), 800, 600);
    hist_m->DrawCopy(); 

    c = new TCanvas("c_dt", Form("data: %s",path_infile.data()), 800, 600); 
    hist_dt->DrawCopy(); 

    return 0; 
}

#endif 