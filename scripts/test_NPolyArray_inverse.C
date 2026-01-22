
//APEX headers
#include "include/RDFNodeAccumulator.h"
#include "include/Get_TParameter_from_TFile.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include <ApexOptics.h> 
//ROOT headers
#include <TParameter.h>
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h> 
#include <TLine.h> 
//stdlib headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>


using namespace std; 
using namespace ROOT::VecOps; 

using ApexOptics::Trajectory_t; 

//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int test_NPolyArray_inverse(    const char* path_infile="data/replay/real_L_V2_sieve.root",  
                                const char* tree_name="tracks_fp" ) 
{
    auto infile = new TFile(path_infile, "READ");

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(__func__, "libApexOptics could not be loaded."); 
        return 1; 
    }
    
    //check if the infile could be opened
    if (!infile || infile->IsZombie()) {
        Error(__func__, "root file '%s' could not be opened.", path_infile); 
        return 1; 
    }
    infile->Close(); 
    delete infile; 

    const bool is_RHRS = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS").value_or(false); 
       
    
    vector<string> branches_sv{
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    }; 

    vector<string> branches_fp{
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 


    cout << "Parsing arrays..." << flush; 

    NPolyArray parr_fp_sv = ApexOptics::Parse_NPolyArray_from_file("data/poly/mc_fp_sv_L_5ord.dat", branches_sv, 4); 
    NPolyArray parr_sv_fp = ApexOptics::Parse_NPolyArray_from_file("data/poly/mc_sv_fp_L_5ord.dat", branches_fp, 5); 
    
    cout << "done." << endl; 
    //now, we're ready to deal with the data. 

    ROOT::EnableImplicitMT(); 
    Info(__func__, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //this is a little helper class which is meant to avoid some of the awkward syntax typically associated with RDataFrames creation. 
    auto rna = RDFNodeAccumulator(df); 

    //plot of x_fp vs reco_x_fp
    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{ x, y, dxdz, dydz };
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Define("Xsv_reco",  [&parr_fp_sv](const Trajectory_t& Xfp)
        {
            return ApexOptics::RVec_to_Trajectory_t(
                parr_fp_sv.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz})
            ); 
        }, {"Xfp"});

    rna.Define("Xfp_reco", [&parr_sv_fp](double x, double y, double dxdz, double dydz, double dpp)
        {
            return ApexOptics::RVec_to_Trajectory_t(
                parr_sv_fp.Eval({x,y,dxdz,dydz,dpp})
            ); 
        }, {"x_sv", "y_sv", "dxdz_sv", "dydz_sv", "dpp_sv"});

    rna.Define("Xfp_error", [](const Trajectory_t& Xfp, const Trajectory_t& Xfp_reco)
        {
            return Xfp + (Xfp_reco*-1.);
        }, {"Xfp", "Xfp_reco"}); 
    

    
    //Error of reconsturctions 

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(__func__, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xfp_reco", {
        {"reco_x_fp",    &Trajectory_t::x},
        {"reco_y_fp",    &Trajectory_t::y},
        {"reco_dxdz_fp", &Trajectory_t::dxdz},
        {"reco_dydz_fp", &Trajectory_t::dydz}
    }); 

    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xfp_error", {
        {"err_x_fp",    &Trajectory_t::x},
        {"err_y_fp",    &Trajectory_t::y},
        {"err_dxdz_fp", &Trajectory_t::dxdz},
        {"err_dydz_fp", &Trajectory_t::dydz}
    }); 

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(__func__, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    //drawing ranges for different coordinates
    constexpr double x_range[]    = {-0.8, +0.8};
    constexpr double y_range[]    = {-0.065, +0.050};
    constexpr double dxdz_range[] = {-0.03, +0.03};
    constexpr double dydz_range[] = {-0.05, +0.05};


    //create both histograms
    auto hist_x = rna.Get()
        .Histo2D<double>({"h_x", "x-fp reco;x_{fp};reco x_{fp}", 
            200, x_range[0],x_range[1], 
            200, x_range[0],x_range[1]}, 
        "x_fp", "reco_x_fp");
    
    //create both histograms
    auto hist_y = rna.Get()
        .Histo2D<double>({"h_y", "y-fp reco;y_{fp};reco y_{fp}", 
            200, y_range[0],y_range[1], 
            200, y_range[0],y_range[1]}, 
        "y_fp", "reco_y_fp");

    //create both histograms
    auto hist_dxdz = rna.Get()
        .Histo2D<double>({"h_dxdz", "dx/dz-fp reco;dx/dz_{fp};reco dx/dz_{fp}", 
            200, dxdz_range[0],dxdz_range[1], 
            200, dxdz_range[0],dxdz_range[1]}, 
        "dxdz_fp", "reco_dxdz_fp");

    //create both histograms
    auto hist_dydz = rna.Get()
        .Histo2D<double>({"h_dydz", "dy/dz-fp reco;dydz_{fp};reco dydz_{fp}", 
            200, dydz_range[0],dydz_range[1], 
            200, dydz_range[0],dydz_range[1]}, 
        "dydz_fp", "reco_dydz_fp");

    
    constexpr double error_x = 5e-3; 
    constexpr double error_y = 5e-3; 
    constexpr double error_dxdz = 5e-3; 
    constexpr double error_dydz = 5e-3; 

    //create both histograms
    auto hist_x_err = rna.Get()
        .Histo1D<double>({"herr_x", "x-fp reco error;x_{fp} - reco x_{fp}", 200, -error_x, +error_x}, "err_x_fp");
    
    //create both histograms
    auto hist_y_err = rna.Get()
        .Histo1D<double>({"herr_y", "y-fp reco error;y_{fp} - reco y_{fp}", 200, -error_y, +error_y}, "err_y_fp");
    
    //create both histograms
    auto hist_dxdz_err = rna.Get()
        .Histo1D<double>({"herr_dxdz", "dxdz-fp reco error;dxdz_{fp} - reco dxdz_{fp}", 200, -error_dxdz, +error_dxdz}, "err_dxdz_fp");
    
    //create both histograms
    auto hist_dydz_err = rna.Get()
        .Histo1D<double>({"herr_dydz", "dydz-fp reco error;dydz_{fp} - reco dydz_{fp}", 200, -error_dydz, +error_dydz}, "err_dydz_fp");
    

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(1); 

    //PolynomialCut::InteractiveApp((TH2*)hist_z_y->Clone("hclone"), "col2", kSunset); 
    //return 0; 
    TCanvas* c; 

    c = new TCanvas("c1", Form("2d plots. data: %s",path_infile), 1200, 700); 
    
    c->Divide(2,2, 0.01,0.01); 
    int c_num=1; 

    const char* draw_option = "col2"; 

    auto line = new TLine; 
    line->SetLineColor(kRed);
    line->SetLineWidth(1); 

    (c->cd(c_num++))->SetLeftMargin(0.15); hist_x   ->DrawCopy(draw_option); line->DrawLine(x_range[0],x_range[0],       x_range[1],x_range[1]);
    (c->cd(c_num++))->SetLeftMargin(0.15); hist_y   ->DrawCopy(draw_option); line->DrawLine(y_range[0],y_range[0],       y_range[1],y_range[1]);
    (c->cd(c_num++))->SetLeftMargin(0.15); hist_dxdz->DrawCopy(draw_option); line->DrawLine(dxdz_range[0],dxdz_range[0], dxdz_range[1],dxdz_range[1]);
    (c->cd(c_num++))->SetLeftMargin(0.15); hist_dydz->DrawCopy(draw_option); line->DrawLine(dydz_range[0],dydz_range[0], dydz_range[1],dydz_range[1]);

    
    c = new TCanvas("c2", Form("1d plots. data: %s",path_infile), 1200, 700); 
    
    c->Divide(2,2, 0.01,0.01); 
    c_num=1; 

    (c->cd(c_num++))->SetLeftMargin(0.15); hist_x_err   ->DrawCopy(draw_option); 
    (c->cd(c_num++))->SetLeftMargin(0.15); hist_y_err   ->DrawCopy(draw_option);
    (c->cd(c_num++))->SetLeftMargin(0.15); hist_dxdz_err->DrawCopy(draw_option);
    (c->cd(c_num++))->SetLeftMargin(0.15); hist_dydz_err->DrawCopy(draw_option);



    return 0; 
}