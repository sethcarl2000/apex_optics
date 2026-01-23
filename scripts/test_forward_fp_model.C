
//APEX headers
#include "include/RDFNodeAccumulator.h"
#include "include/TestAngleReco.h"
#include "include/ChainedOpticsModel.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include "include/Get_TParameter_from_TFile.h"
#include <NPolyArrayChain.h> 
#include <Interactive3dHist.hxx>
#include <ApexOptics.h> 
//ROOT headers
#include <ApexOptics.h> 
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h> 
#include <TGClient.h> 
#include <ROOT/RVec.hxx>
//stdlib headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <vector> 
#include <map>
#include <utility> 
#include <limits> 

using namespace std; 
using namespace ROOT::VecOps; 
using ApexOptics::OpticsTarget_t; 
using ApexOptics::Trajectory_t; 

namespace {
        
    constexpr double range_dxdz[] = {-0.05, +0.05}; 
    constexpr double range_dydz[] = {-0.04, +0.04}; 

    const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    const vector<string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"};

    const vector<string> branches_fwd_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
    const vector<string> branches_fwd_fp{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"};

    const vector<string> branches_rev_sv{"fwd_x_sv","fwd_y_sv","fwd_dxdz_sv","fwd_dydz_sv","fwd_dpp_sv"};
    const vector<string> branches_rev_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
}; 

//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int test_forward_fp_model( const char* path_infile ="data/replay/real_L_V2.root",
                           const char* target_name ="V2" ) 
{
    const char* tree_name = "tracks_fp";

    auto infile = new TFile(path_infile, "READ");

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(__func__, "libApexOptics could not be loaded."); 
        return -1; 
    }
    
    const bool is_RHRS = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS").value_or(false);    

    //try to get the target we need. 
    OpticsTarget_t target; 
    try { 
        target = ApexOptics::GetTarget(string(target_name)); 
    
    } catch (const std::exception& e) {

        Error(__func__, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
        return -1; 
    }

    //in the syntax below, an '<=' arrow represents an input or output of a polynomial. 
    // if a polynomial is written '[Poly]' then it is trained soley on monte-carlo data. 
    // if a polynomial is written '_Poly_', then it is trained on real data, which is fed into monte-carlo polynomials 
    //      (in order to reconstruct XQ1).
    // uncomment whichever configuration you want to use.
    // 

    NPolyArray parr_fp_fpfwd = ApexOptics::Parse_NPolyArray_from_file(
        "data/poly/mc_fp_fp-fwd_L_5ord.dat", 
        {"fwd_x_fp", "fwd_y_fp", "fwd_dxdz_fp", "fwd_dydz_fp"}, 
        4
    ); 

    ChainedOpticsModel* model = new ChainedOpticsModel(is_RHRS); 
    model->CreateChainRev({

        // sv <= _Poly_ <= fp
        {"data/poly/mc_fp_sv_L_5ord.dat", branches_sv, 4} 
    }); 

    model->CreateChainFwd({
    
        // sv => _Poly_ => fp 
        {"data/poly/mc_sv_fp_L_5ord.dat", branches_fp, 5},
    }); 

    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 

    ROOT::EnableImplicitMT(); 
    Info(__func__, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //this is a little helper class which is meant to avoid some of the awkward syntax typically associated with RDataFrames creation. 
    auto rna = RDFNodeAccumulator(tree_name, path_infile); 

    //probably not the most elegant way to do this, but __func__ we are. 
    rna.Define("Xfp_fwd", [&parr_fp_fpfwd](double x, double y, double dxdz, double dydz)
        {
            return ApexOptics::RVec_to_Trajectory_t(parr_fp_fpfwd.Eval({x,y,dxdz,dydz})); 
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Define("Xsv_first_guess", [&model](Trajectory_t Xfp_fwd)
        {
            return model->Compute_Xsv_first_guess(Xfp_fwd);
        }, {"Xfp_fwd"});

    
    const int n_iterations = 3;  
    const double fp_error_threshold = 1e-6; 

    rna.Define("Xsv_reco", [&model, n_iterations, fp_error_threshold](Trajectory_t Xfp_fwd, Trajectory_t Xsv, TVector3 vtx_hcs)
        {
            //yes, i know this syntax is disgusting. sorry. I will (hopefully) fix it sometime soon. 
            NPolyArray* arr_sv_fp = &model->GetChain_fwd()->arrays[0].first; 

            auto&& Xsv_rvec = ApexOptics::Trajectory_t_to_RVec(Xsv); 
            auto&& Xfp_rvec = ApexOptics::Trajectory_t_to_RVec(Xfp_fwd); 
            
            for (int n_it=0; n_it<n_iterations; n_it++) {
                
                //first, use the first-order taylor expansion of our optics model to match-up with the vertical raster
                Xsv = model->Compute_Xsv(Xfp_fwd, vtx_hcs, Xsv_rvec);

                //now, match up the new, interpolated value of Xsv to the known value of Xfp.  
                Xsv_rvec = ApexOptics::Trajectory_t_to_RVec(Xsv); 
                arr_sv_fp->Iterate_to_root(Xsv_rvec, Xfp_rvec, 10, fp_error_threshold); 
                Xsv = ApexOptics::RVec_to_Trajectory_t(Xsv_rvec); 
            }
            return Xsv;

        }, {"Xfp_fwd", "Xsv_first_guess", "position_vtx"});
        
    rna.Define("Xhcs_reco", [is_RHRS](const Trajectory_t& Xsv)
        {
            return ApexOptics::SCS_to_HCS(is_RHRS, Xsv); 
        }, {"Xsv_reco"});

    rna.Overwrite("z_reco_horizontal", [](const Trajectory_t& Xhcs, const TVector3& vtx)
        {   
            return - ( Xhcs.y - vtx.y() ) / Xhcs.dydz; 
        }, {"Xhcs_reco", "position_vtx"});

    rna.Overwrite("z_reco_vertical",   [](const Trajectory_t& Xhcs, const TVector3& vtx)
        {
            return - ( Xhcs.x - vtx.x() ) / Xhcs.dxdz; 
        }, {"Xhcs_reco", "position_vtx"});

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(__func__, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }
    
    rna = Add_branches_from_Trajectory_t(rna.Get(), "Xsv_first_guess", {
        "fg_x_sv",
        "fg_y_sv",
        "fg_dxdz_sv",
        "fg_dydz_sv"
    }); 

    rna = Add_branches_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        "reco_x_sv",
        "reco_y_sv",
        "reco_dxdz_sv",
        "reco_dydz_sv"
    }); 

    rna.Define("y_pos", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}); 

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(__func__, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    const double min_cerenkov_sum = 1.5e3; 
    const double min_Esh_Eps_sum = 0.5; 

    
    auto hist_z_dy_nocut = rna.Get()
        .Histo2D<double>({"h_z_dy_nocut", "No PID cut;z_{tg};dy/dz_{sv}", 200, -0.4, +0.4, 200, -0.04, 0.03}, "z_reco_vertical", "reco_dydz_sv"); 

   

    auto hist_z = rna.Get()
        .Histo1D<double>({"h_z", "Z_{hcs} reconstruction;z_{tg};", 200, -0.4, +0.4}, "z_reco_vertical"); 

    
    //create both histograms
    auto hist_xy = rna.Get()
        .Histo2D({"h_xy",     "x_{sv} vs y_{sv} (linear Jacobian-correction);x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "reco_x_sv", "reco_y_sv");
    
    auto hist_angles = rna.Get()
        .Histo2D({"h_angles", "dx/dz_{sv} vs dy/dz_{sv} (linear Jacobian-correction);dx/dx_{sv};dy/dz_{sv}", 200, range_dxdz[0],range_dxdz[1], 200, range_dydz[0],range_dydz[1]}, "reco_dxdz_sv", "reco_dydz_sv"); 
    

    //create both histograms
    auto hist_fg_xy = rna.Get()
        .Histo2D({"h_fg_xy",     "x_{sv} vs y_{sv};x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "fg_x_sv", "fg_y_sv");

    auto hist_fg_angles = rna.Get()
        .Histo2D({"h_fg_angles", "dx/dz_{sv} vs dy/dz_{sv};dx/dx_{sv};dy/dz_{sv}", 200, range_dxdz[0],range_dxdz[1], 200, range_dydz[0],range_dydz[1]}, "fg_dxdz_sv", "fg_dydz_sv");

    char c_title[255]; 
    sprintf(c_title, "data:'%s' target: %s", path_infile, target.name.c_str()); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 

    /*/Measure lab-vertical angle (dxdz)
    TestAngleReco::Evaluate(is_RHRS, TestAngleReco::kDxdz, rna.Get(), target, 
        4,13, 
        1,8, 
        1.35, 0.50, 1.00,
        "reco_dxdz_sv", "reco_dydz_sv",
        false, 
        range_dxdz[0], range_dxdz[1],
        range_dydz[0], range_dydz[1]
    ); 

    //Measure lab-horizontal angle (dydz)
    TestAngleReco::Evaluate(is_RHRS, TestAngleReco::kDydz, rna.Get(), target, 
        4,12, 
        1,10, 
        1.25, 0.65, 0.20,
        "reco_dxdz_sv", "reco_dydz_sv",
        true, 
        range_dxdz[0], range_dxdz[1],
        range_dydz[0], range_dydz[1]
    ); 
    return 0;  //*/ 

    TCanvas *c; 
    c = new TCanvas("c_z_reco", c_title, 700, 700); 
    hist_z->DrawCopy(); 


    c = new TCanvas("c1", c_title, 1200, 600); 
    c->SetLeftMargin(0.12); c->SetRightMargin(0.05); 
    c->Divide(2,1, 0.01,0.01); 

    c->cd(1)->SetLeftMargin(0.15); hist_fg_xy->DrawCopy("col2"); 
    c->cd(2)->SetLeftMargin(0.15); hist_xy->DrawCopy("col2");
    
    c = new TCanvas("c2", c_title, 1200, 600); 
    c->Divide(2,1, 0.01,0.01); 

    c->cd(1)->SetLeftMargin(0.15); hist_fg_angles->DrawCopy("col2"); 
    c->cd(2)->SetLeftMargin(0.15); hist_angles->DrawCopy("col2"); 

    return 0; 
}