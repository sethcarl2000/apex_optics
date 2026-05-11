#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <vector> 
#include <map>
#include <utility> 
#include "include/RDFNodeAccumulator.h"
#include "include/TestAngleReco.h"
#include "include/ChainedOpticsModel.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include <TParameter.h>
#include <ApexOptics.h> 
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h> 
#include <NPolyArrayChain.h> 
#include <TRandom3.h> 
#include <limits> 
#include <Interactive3dHist.hxx>
#include <TGClient.h> 

using namespace std; 
using namespace ROOT::VecOps; 
using ApexOptics::OpticsTarget_t; 
using ApexOptics::Trajectory_t; 

namespace {
    constexpr double pi = 3.14159265359; 
    
    //z-heights of different VDC planes (in meters)
    constexpr double z_U1 = 0.;     
    constexpr double z_V1 = 26e-3; 
    constexpr double z_U2 = 0.335; 
    constexpr double z_V2 = 0.335 + 26e-3; 

    TRandom3 global_rand; 
}

Trajectory_t Add_VDC_smearing(Trajectory_t Xfp, double VDC_smearing)  
{
    //convert to 'transport' coordinates, in which only 'theta' (dx/dz) is different
    Trajectory_t Xtr = Xfp; 
    Xtr.dxdz += Xtr.x/6.; 

    //now, rotate to 'det' coordinates
    TVector3 s( Xtr.dxdz, Xtr.dydz, 1. ); //'direction' vector
    TVector3 r( Xtr.x,    Xtr.y,    0. ); //'position'  vector 

    s.RotateY( pi/4. ); s.RotateZ( pi/4. ); 
    r.RotateY( pi/4. ); r.RotateZ( pi/4. ); 

    double ru = r[0]; 
    double rv = r[1]; 
    double rz = r[2]; 

    double su = s[0]/s[2]; 
    double sv = s[1]/s[2]; 

    //compute the 4 VDC intercepts 
    double u1 = ru + (z_U1 - rz)*su + global_rand.Gaus()*VDC_smearing; 
    double v1 = rv + (z_V1 - rz)*sv + global_rand.Gaus()*VDC_smearing; 
    double u2 = ru + (z_U2 - rz)*su + global_rand.Gaus()*VDC_smearing; 
    double v2 = rv + (z_V2 - rz)*sv + global_rand.Gaus()*VDC_smearing;

    //now, re-compute the angle and displacement TVector3's components
    su = (u2 - u1)/(z_U2 - z_U1); 
    sv = (v2 - v1)/(z_V2 - z_V1); 

    ru = u1 - (z_U1 - rz)*su; 
    rv = v1 - (z_V1 - rz)*sv; 

    //now, de-rotate the new vector (with smearing added).
    s = TVector3( su, sv, 1. ); 
    r = TVector3( ru, rv, rz ); 

    s.RotateZ( -pi/4. ); s.RotateY( -pi/4. ); 
    r.RotateZ( -pi/4. ); r.RotateY( -pi/4. ); 
    
    //now, we have our vectors back in 'transport' coordinates. let's convert to a 'Trajectory_t' object: 
    Xtr = Trajectory_t{
        .x = r.x() + ( 0.-r.z() )*(s.x()/s.z()), 
        .y = r.y() + ( 0.-r.z() )*(s.y()/s.z()), 
        .dxdz = s.x()/s.z(),
        .dydz = s.y()/s.z()
    }; 

    //now, perform the conversion from 'transport' to 'focal-plane' coordinates
    Xfp = Xtr; 
    Xfp.dxdz += -Xfp.x/6.; 

    return Xfp; 
}


const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
const vector<string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"};

const vector<string> branches_fwd_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
const vector<string> branches_fwd_fp{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"};

const vector<string> branches_rev_sv{"fwd_x_sv","fwd_y_sv","fwd_dxdz_sv","fwd_dydz_sv","fwd_dpp_sv"};
const vector<string> branches_rev_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};

//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int simulate_VDC_smearing(  const char* path_infile ="data/replay/real_L_V1-dp.root",
                            const char* target_name ="V1",
                            const double VDC_smearing = 85e-6, 
                            const char* path_dbfile ="data/csv/poly_WireAndFoil_fp_sv_L_4ord.dat",  
                            const char* tree_name   ="tracks_fp" ) 
{

    const char* const here = "test_forward_chain"; 

    auto infile = new TFile(path_infile, "READ");

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(here, "libApexOptics could not be loaded."); 
        return 1; 
    }
    
    //check if the infile could be opened
    if (!infile || infile->IsZombie()) {
        Error(here, "root file '%s' could not be opened.", path_infile); 
        return 1; 
    }

    //check if we can find the proper TTree
    TTree* tree = (TTree*)infile->Get(tree_name); 
    if (!tree) {
        Error(here, "could not find TTree '%s'", tree_name); 
        return 1; 
    }

    //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
    TParameter<bool>* param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
    if (!param_is_RHRS) {
        Error(here, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
        return 1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 
       
    delete tree; 
    infile->Close(); 
    delete infile; 

    //try to get the target we need. 
    OpticsTarget_t target; 
    try { 
        target = ApexOptics::GetTarget(string(target_name)); 
    
    } catch (const std::exception& e) {

        Error(here, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
        return -1; 
    }

    //in the syntax below, an '<=' arrow represents an input or output of a polynomial. 
    // if a polynomial is written '[Poly]' then it is trained soley on monte-carlo data. 
    // if a polynomial is written '_Poly_', then it is trained on real data, which is fed into monte-carlo polynomials 
    //      (in order to reconstruct XQ1).
    // uncomment whichever configuration you want to use. 
    // 

    struct NPolyArrayConstructor_t { string path{}; vector<string> coords{}; int input_DoF{0}; }; 


    ChainedOpticsModel* model = new ChainedOpticsModel(is_RHRS); 
    model->CreateChainRev({

        // sv <= [Poly] <= fp
        //{"data/poly/fits_5Dec/V123_fp_sv_4ord.dat", branches_sv, 4} //*/ 
        {"data/poly/fits_21Dec/V123_fp_sv_L_4ord.dat", branches_sv, 4} //*/  
    }); 
    model->CreateChainFwd({
        // sv => [Poly] => fp-fwd => _Poly_ => fp 
        {"data/csv/poly_prod_sv_fp_L_4ord.dat",         branches_fp, 5},
        {"data/poly/fits_6Nov/V123_fp-fwd_fp_1ord.dat", branches_fp, 4} 
         //*/ 
    
    }); 

    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame *df = nullptr;    
    try { 

        df = new ROOT::RDataFrame(tree_name, path_infile);

    } catch (const std::exception& e) {

        Error (here, "Something went wrong trying to create the RDataFrame.\n what(): %s", e.what()); 
        return -1; 
    } 

    //this is a little helper class which is meant to avoid some of the awkward syntax typically associated with RDataFrames creation. 
    auto rna = RDFNodeAccumulator(*df); 

    //probably not the most elegant way to do this, but here we are. 
    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{ x, y, dxdz, dydz };
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Overwrite("Xsv_first_guess", [&model](Trajectory_t Xfp)
        {
            return model->Compute_Xsv_first_guess(Xfp); 
        }, {"Xfp"});

    rna.Overwrite("Xsv_reco", [&model](Trajectory_t Xfp, TVector3 vtx_hcs)
        {
            return model->Compute_Xsv(Xfp, vtx_hcs);
        }, {"Xfp", "position_vtx"});

    rna.Define("Xfp_smeared", [VDC_smearing](Trajectory_t Xfp)
        {
            Trajectory_t Xfp_smeared = Add_VDC_smearing(Xfp, VDC_smearing); 

            return Xfp_smeared; 
        }, {"Xfp"});

    rna.Define("Xfp_smearing_error", [](Trajectory_t Xfp_smeared, Trajectory_t Xfp)
        {
            return (Xfp_smeared - Xfp)*1e3; 
        }, {"Xfp_smeared", "Xfp"}); 

    //now, compute the 'smeared' version of the Xsv reconstructed target coordinate vector
    rna.Overwrite("Xsv_reco_smeared", [&model](Trajectory_t Xfp_smeared, TVector3 vtx_hcs)
        {
            return model->Compute_Xsv(Xfp_smeared, vtx_hcs);
        }, {"Xfp_smeared", "position_vtx"});
    rna.Overwrite("Xsv_smearing_error", [&model](Trajectory_t Xsv, Trajectory_t Xsv_smeared)
        {
            return (Xsv_smeared - Xsv)*1e3; 
        }, {"Xsv_reco", "Xsv_reco_smeared"});

    
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
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }
    
    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xsv_first_guess", {
        {"fg_x_sv",    &Trajectory_t::x},
        {"fg_y_sv",    &Trajectory_t::y},
        {"fg_dxdz_sv", &Trajectory_t::dxdz},
        {"fg_dydz_sv", &Trajectory_t::dydz}
    }); 

    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        {"reco_x_sv",    &Trajectory_t::x},
        {"reco_y_sv",    &Trajectory_t::y},
        {"reco_dxdz_sv", &Trajectory_t::dxdz},
        {"reco_dydz_sv", &Trajectory_t::dydz}
    }); 

    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xfp_smearing_error", {
        {"err_x_fp",    &Trajectory_t::x},
        {"err_y_fp",    &Trajectory_t::y},
        {"err_dxdz_fp", &Trajectory_t::dxdz},
        {"err_dydz_fp", &Trajectory_t::dydz}
    }); 
    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xsv_smearing_error", {
        {"err_x_sv",    &Trajectory_t::x},
        {"err_y_sv",    &Trajectory_t::y},
        {"err_dxdz_sv", &Trajectory_t::dxdz},
        {"err_dydz_sv", &Trajectory_t::dydz},
        {"err_dpp_sv",  &Trajectory_t::dpp}
    }); 

    rna.Define("y_pos", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}); 

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    //range of fp-drawing in mrad 
    const double err_bound_fp = 2.5; 

    auto hist_x_fp = rna.Get() 
        .Histo1D<double>({"h_xfp", "VDC-smearing error of x_{fp};error (mm)", 200, -err_bound_fp, +err_bound_fp}, "err_x_fp"); 
    auto hist_y_fp = rna.Get() 
        .Histo1D<double>({"h_yfp", "VDC-smearing error of y_{fp};error (mm)", 200, -err_bound_fp, +err_bound_fp}, "err_y_fp"); 
    auto hist_dxdz_fp = rna.Get() 
        .Histo1D<double>({"h_dxdzfp", "VDC-smearing error of dx/dz_{fp};error (mrad)", 200, -err_bound_fp, +err_bound_fp}, "err_dxdz_fp"); 
    auto hist_dydz_fp = rna.Get() 
        .Histo1D<double>({"h_dydzfp", "VDC-smearing error of dy/dz_{fp};error (mrad)", 200, -err_bound_fp, +err_bound_fp}, "err_dydz_fp"); 

    auto hist_x_sv = rna.Get() 
        .Histo1D<double>({"h_xfp", "VDC-smearing error of x_{sv};error (mm)", 200, -err_bound_fp, +err_bound_fp}, "err_x_sv"); 
    auto hist_y_sv = rna.Get() 
        .Histo1D<double>({"h_yfp", "VDC-smearing error of y_{sv};error (mm)", 200, -err_bound_fp, +err_bound_fp}, "err_y_sv"); 
    auto hist_dxdz_sv = rna.Get() 
        .Histo1D<double>({"h_dxdzfp", "VDC-smearing error of dx/dz_{sv};error (mrad)", 200, -err_bound_fp, +err_bound_fp}, "err_dxdz_sv"); 
    auto hist_dydz_sv = rna.Get() 
        .Histo1D<double>({"h_dydzfp", "VDC-smearing error of dy/dz_{sv};error (mrad)", 200, -err_bound_fp, +err_bound_fp}, "err_dydz_sv"); 


    //create both histograms
    auto hist_xy = rna.Get()
        .Histo2D({"h_xy",     "x_{sv} vs y_{sv} (linear dp/p correction);x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "reco_x_sv", "reco_y_sv");
    
    auto hist_angles = rna.Get()
        .Histo2D({"h_angles", "dx/dz_{sv} vs dy/dz_{sv} (linear dp/p correction);dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "reco_dxdz_sv", "reco_dydz_sv"); 
    

    //create both histograms
    auto hist_fg_xy = rna.Get()
        .Histo2D({"h_fg_xy",     "x_{sv} vs y_{sv};x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "fg_x_sv", "fg_y_sv");

    auto hist_fg_angles = rna.Get()
        .Histo2D({"h_fg_angles", "dx/dz_{sv} vs dy/dz_{sv};dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "fg_dxdz_sv", "fg_dydz_sv");

    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    //gStyle->SetOptStat(0); 


    TCanvas *c; 
    c = new TCanvas("c_fp_err", c_title, 1600, 900); 
    gPad->SetLeftMargin(0.15); 
    c->Divide(2,2);
    c->cd(1); hist_x_fp->DrawCopy(); 
    c->cd(2); hist_y_fp->DrawCopy(); 
    c->cd(3); hist_dxdz_fp->DrawCopy(); 
    c->cd(4); hist_dydz_fp->DrawCopy(); 

    c = new TCanvas("c_sv_err", c_title, 1600, 900); 
    gPad->SetLeftMargin(0.15); 
    c->Divide(2,2);
    c->cd(1); hist_x_sv->DrawCopy(); 
    c->cd(2); hist_y_sv->DrawCopy(); 
    c->cd(3); hist_dxdz_sv->DrawCopy(); 
    c->cd(4); hist_dydz_sv->DrawCopy(); 

    return 0; 
}