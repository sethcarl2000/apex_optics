
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include "include/RDFNodeAccumulator.h"
#include <TParameter.h>
#include <ApexOptics.h> 
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h> 
#include <TRandom3.h> 


using namespace std; 
using namespace ROOT::VecOps; 

using ApexOptics::Trajectory_t; 

//_______________________________________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Trajectory_t struct. 
ROOT::RDF::RNode add_branch_from_Trajectory_t(ROOT::RDF::RNode df, const char* branch_in, map<string, double Trajectory_t::*> branches)
{
    const int n_nodes = branches.size() + 1; 
    RVec<ROOT::RDF::RNode> nodes{ df }; 

    int i_branch=0; 
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double Trajectory_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Trajectory_t::coord' 
        nodes.push_back( nodes.back()
            
            .Define(branch_name, [coord](const Trajectory_t& track) { return track.*coord; }, {branch_in})
        ); 
    }
    
    return nodes.back(); 
}


Trajectory_t RVec_to_Trajectory_t(const RVec<double>& v){
    
    //check if input vector is right size. 
    if (!(v.size() == 4 || v.size() == 5)) {
        ostringstream oss; 
        oss << "in <RVec_to_Trajectory_t>: vector given is wrong size (" << v.size() << "), "
               " must be either 4 or 5."; 
        throw invalid_argument(oss.str()); 
        return Trajectory_t{}; 
    }

    return Trajectory_t{
        .x=v[0], 
        .y=v[1], 
        .dxdz=v[2],
        .dydz=v[3],
        .dpp= (v.size() == 4 ? std::numeric_limits<double>::quiet_NaN() : v[4])
    }; 
}


//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int test_linear_model( const char* path_infile="data/replay/real_L_V2.root",
                        const char* path_dbfile="data/csv/poly_WireAndFoil_L_4ord.dat",  
                        const char* path_dXsv="data/csv/poly_dXsv_L_3ord.dat",
                        const char* tree_name="tracks_fp" ) 
{
    const char* const here = "test_linear_model"; 

    //get the arm we need. 
    TFile *infile; 
    TParameter<bool> *param_is_RHRS; 

    try {
        infile = new TFile(path_infile, "READ");

        //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
        param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
        if (!param_is_RHRS) {
            Error(here, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
            return 1; 
        }
    
    } catch (const std::exception& e) {

        Error(here, "Execption caught trying to open file.\n what(): %s", e.what()); 
        return -1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 
       
    infile->Close(); 
    delete infile; 

    
    const vector<string> branches_sv{
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    }; 

    const vector<string> branches_fp{
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    const vector<string> branches_dXsv{
        "d_x_sv",
        "d_y_sv",
        "d_dxdz_sv",
        "d_dydz_sv"
    }; 

    NPolyArray parr_fp_sv, parr_dXsv;  

    try { 
        parr_fp_sv = ApexOptics::Parse_NPolyArray_from_file(path_dbfile, branches_sv, 4); 
        parr_dXsv  = ApexOptics::Parse_NPolyArray_from_file(path_dXsv, branches_dXsv, 5); 
    
    } catch (const std::exception& e) {

        Error(here, "Error caught trying to parse NPolyArray from file.\n what(): %s", e.what()); 
        return -1; 
    }

    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 

    const double y_rast_rms = 0.113e-3; 

    const double d_dpp = 1e-3; 

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //this is a little helper class which is meant to avoid some of the awkward syntax typically associated with RDataFrames creation. 
    auto rna = RDFNodeAccumulator(df); 

    //this branch would only be defined if this is monte-carlo data. 
    const bool is_MonteCarlo = false; //rna.IsBranchDefined("x_sv"); 

    cout << (is_MonteCarlo ? "Monte-carlo mode" : "Real-data mode") << endl; 

    TRandom3 rand; 

    //probably not the most elegant way to do this, but here we are. 
    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{ .x=x, .y=y, .dxdz=dxdz, .dydz=dydz }; 

        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Define("Xsv_first_guess",  [&parr_fp_sv](const Trajectory_t& Xfp)
        {
            return RVec_to_Trajectory_t( parr_fp_sv.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz}) ); 

        }, {"Xfp"});

    //define a y-react pos. if it's monte-carlo data, then we want to add a realistic RMS. 
    rna.Define("y_vtx", [y_rast_rms, &rand, is_MonteCarlo](TVector3 vtx)
        {
            return vtx.y() + ( is_MonteCarlo ? rand.Gaus() * y_rast_rms : 0. ); 

        }, {"position_vtx"}); 

    rna.Define("dXsv",     [&parr_dXsv](Trajectory_t Xsv)
        {
            return RVec_to_Trajectory_t( parr_dXsv.Eval({Xsv.x, Xsv.y, Xsv.dxdz, Xsv.dydz, Xsv.dpp}) ); 
        }, {"Xsv_first_guess"}); 
    
    auto Compute_dp = [d_dpp, is_RHRS](Trajectory_t Xsv, Trajectory_t dXsv, double y_vtx, double z_vtx)
    {
        Trajectory_t Xsv_dp{
                Xsv.x       + d_dpp * dXsv.x,
                Xsv.y       + d_dpp * dXsv.y, 
                Xsv.dxdz    + d_dpp * dXsv.dxdz, 
                Xsv.dydz    + d_dpp * dXsv.dydz, 
                Xsv.dpp     + d_dpp
            }; 

            Trajectory_t Xhcs    = ApexOptics::SCS_to_HCS(is_RHRS, Xsv); 
            Trajectory_t Xhcs_dp = ApexOptics::SCS_to_HCS(is_RHRS, Xsv_dp); 

            // now, project the track onto the target's z-plane.
            double y1 = Xhcs_dp.y  +  Xhcs_dp.dydz * z_vtx; 
            double y0 = Xhcs.y     +  Xhcs.dydz    * z_vtx;

            // y_vtx = y0 + sy * dp;        sy := dy/dp = ( y1 - y0 )/( p1 - p0 ); 

            return ( y_vtx - y0 ) * d_dpp / ( y1 - y0 ); 
    };


    rna.Define("dp", [&Compute_dp](Trajectory_t Xsv, Trajectory_t dXsv, double y_vtx, TVector3 vtx)
        {   
            return Compute_dp(Xsv, dXsv, y_vtx, vtx.z()); 

        }, {"Xsv_first_guess", "dXsv", "y_vtx", "position_vtx"});

    rna.Define("dp_upfoil", [&Compute_dp](Trajectory_t Xsv, Trajectory_t dXsv, double y_vtx, TVector3 vtx)
        {   
            return Compute_dp(Xsv, dXsv, y_vtx, vtx.z() + 55e-3); 

        }, {"Xsv_first_guess", "dXsv", "y_vtx", "position_vtx"});

    rna.Define("dp_downfoil", [&Compute_dp](Trajectory_t Xsv, Trajectory_t dXsv, double y_vtx, TVector3 vtx)
        {   
            return Compute_dp(Xsv, dXsv, y_vtx, vtx.z() - 55e-3); 

        }, {"Xsv_first_guess", "dXsv", "y_vtx", "position_vtx"});
 

    rna.Define("Xsv_reco", [](Trajectory_t Xsv, Trajectory_t dXsv, double dp)
        {
            return Trajectory_t{
                Xsv.x       + dp * dXsv.x,
                Xsv.y       + dp * dXsv.y, 
                Xsv.dxdz    + dp * dXsv.dxdz, 
                Xsv.dydz    + dp * dXsv.dydz, 
                Xsv.dpp     + dp
            }; 

        }, {"Xsv_first_guess", "dXsv", "dp"});

    rna.Define("Xsv_reco_upfoil", [](Trajectory_t Xsv, Trajectory_t dXsv, double dp)
        {
            return Trajectory_t{
                Xsv.x       + dp * dXsv.x,
                Xsv.y       + dp * dXsv.y, 
                Xsv.dxdz    + dp * dXsv.dxdz, 
                Xsv.dydz    + dp * dXsv.dydz, 
                Xsv.dpp     + dp
            }; 

        }, {"Xsv_first_guess", "dXsv", "dp_upfoil"});

    rna.Define("Xsv_reco_downfoil", [](Trajectory_t Xsv, Trajectory_t dXsv, double dp)
        {
            return Trajectory_t{
                Xsv.x       + dp * dXsv.x,
                Xsv.y       + dp * dXsv.y, 
                Xsv.dxdz    + dp * dXsv.dxdz, 
                Xsv.dydz    + dp * dXsv.dydz, 
                Xsv.dpp     + dp
            }; 

        }, {"Xsv_first_guess", "dXsv", "dp_downfoil"});

    rna.Define("Xhcs_reco", [is_RHRS](const Trajectory_t& Xsv)
        {
            return ApexOptics::SCS_to_HCS(is_RHRS, Xsv);  

        }, {"Xsv_reco"});
    
    rna.Define("Xhcs_reco_upfoil", [is_RHRS](const Trajectory_t& Xsv)
        {
            return ApexOptics::SCS_to_HCS(is_RHRS, Xsv);  

        }, {"Xsv_reco_upfoil"});

    rna.Define("Xhcs_reco_downfoil", [is_RHRS](const Trajectory_t& Xsv)
        {
            return ApexOptics::SCS_to_HCS(is_RHRS, Xsv);  

        }, {"Xsv_reco_downfoil"});

    rna.Overwrite("z_reco_horizontal", [](const Trajectory_t& Xhcs, double y_vtx)
        {   
            return - ( Xhcs.y - y_vtx ) / Xhcs.dydz; 

        }, {"Xhcs_reco", "y_vtx"});

    rna.Overwrite("z_reco_vertical",   [](const Trajectory_t& Xhcs, const TVector3& vtx)
        {
            return - ( Xhcs.x - vtx.x() ) / Xhcs.dxdz; 

        }, {"Xhcs_reco", "position_vtx"});

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    rna = add_branch_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        {"reco_x_sv",    &Trajectory_t::x},
        {"reco_y_sv",    &Trajectory_t::y},
        {"reco_dxdz_sv", &Trajectory_t::dxdz},
        {"reco_dydz_sv", &Trajectory_t::dydz}
    }); 

    rna = add_branch_from_Trajectory_t(rna.Get(), "Xsv_first_guess", {
        {"fg_x_sv",    &Trajectory_t::x},
        {"fg_y_sv",    &Trajectory_t::y},
        {"fg_dxdz_sv", &Trajectory_t::dxdz},
        {"fg_dydz_sv", &Trajectory_t::dydz}
    }); 

    rna.Define("x_hcs_upfoil",   [](Trajectory_t Xhcs){ return Xhcs.x; }, {"Xhcs_upfoil"});
    rna.Define("x_hcs",          [](Trajectory_t Xhcs){ return Xhcs.x; }, {"Xhcs"});
    rna.Define("x_hcs_downfoil", [](Trajectory_t Xhcs){ return Xhcs.x; }, {"Xhcs_downfoil"});
    

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    //create both histograms
    auto hist_xy = rna.Get()
        .Histo2D({"h_xy",     "linear dp/p correction: x_{sv} vs y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "reco_x_sv", "reco_y_sv");
    
    auto hist_angles = rna.Get()
        .Histo2D({"h_angles", "linear dp/p correction: dx/dx_{sv} vs dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "reco_dxdz_sv", "reco_dydz_sv"); 
    

    //create both histograms
    auto hist_fg_xy = rna.Get()
        .Histo2D({"h_fg_xy",     "fp=>sv model: x_{sv} vs y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "fg_x_sv", "fg_y_sv");
    
    auto hist_fg_angles = rna.Get()
        .Histo2D({"h_fg_angles", "fp=>sv model: dx/dx_{sv} vs dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "fg_dxdz_sv", "fg_dydz_sv"); 

    auto hist_dp = rna.Get()
        .Histo1D({"h_dp", "dp/p fix; dp/p - dp/p_{guesss}", 200, -0.003, 0.003}, "dp"); 

    auto hist_dp_upfoil = rna.Get()
        .Histo1D({"h_dp", "dp/p fix: z + 55mm; dp/p - dp/p_{guesss}", 200, -0.003, 0.003}, "dp_upfoil"); 

    auto hist_dp_downfoil = rna.Get()
        .Histo1D({"h_dp", "dp/p fix: z - 55mm; dp/p - dp/p_{guesss}", 200, -0.003, 0.003}, "dp_downfoil"); 

    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 

    //PolynomialCut::InteractiveApp((TH2*)hist_z_y->Clone("hclone"), "col2", kSunset); 
    //return 0; 

    auto c0 = new TCanvas("c0", c_title); 
    hist_dp->DrawCopy(); 

    hist_dp_upfoil->SetFillStyle(3004); 
    hist_dp_upfoil->SetFillColor(kRed);
    hist_dp_upfoil->SetLineColor(kRed);
    hist_dp_upfoil->SetStats(0);  
    hist_dp_upfoil->DrawCopy("SAME");

    hist_dp_downfoil->SetFillStyle(3005); 
    hist_dp_downfoil->SetFillColor(kBlue);
    hist_dp_downfoil->SetLineColor(kBlue);
    hist_dp_downfoil->SetStats(0);  
    hist_dp_downfoil->DrawCopy("SAME");

    gPad->BuildLegend(); 

    return 0; 

    auto c1 = new TCanvas("c1", c_title, 1200, 600); 
    c1->Divide(2,1, 0.01,0.01); 

    c1->cd(1); hist_fg_xy->DrawCopy("col2"); 
    c1->cd(2); hist_xy->DrawCopy("col2");

    auto c2 = new TCanvas("c4", c_title, 1200, 600); 
    c2->Divide(2,1, 0.01,0.01);

    c2->cd(1); hist_fg_angles->DrawCopy("col2");
    c2->cd(2); hist_angles->DrawCopy("col2"); 
    
    return 0; 

    return 0; 
}