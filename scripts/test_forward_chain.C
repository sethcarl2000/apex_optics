
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <vector> 
#include <map>
#include <utility> 
#include "include/RDFNodeAccumulator.h"
#include "include/TestAngleReco.h"
#include <TParameter.h>
#include <ApexOptics.h> 
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h> 
#include <NPolyArrayChain.h> 
#include <limits> 

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
int test_forward_chain( const char* path_infile="data/replay/real_L_V2.root",
                        const char* target_name="V2",
                        const char* path_dbfile="data/csv/poly_WireAndFoil_fp_sv_L_4ord.dat",  
                        const char* tree_name="tracks_fp" ) 
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

    const vector<NPolyArrayConstructor_t> path_and_coords_rev{

        /*/ sv <= [Poly] <= fp
        {"data/csv/poly_prod_fp_sv_L_4ord.dat", branches_sv, 4} //*/ 

        /*/ sv <= [Poly] <= fp-fwd <= _Poly_ <= fp
        {"data/csv/poly_fits_fp_fp-fwd_L_4ord.dat", branches_fwd_fp, 4}, 
        {"data/csv/poly_prod_fp_sv_L_4ord.dat",     branches_sv,     4} //*/ 

        /*/ sv <= [Poly] q1 <= [Poly] <= fp-fwd <= _Poly_ <= fp 
        {"data/csv/poly_fits_fp_fp-fwd_L_4ord.dat", branches_fwd_fp, 4}, 
        {"data/csv/poly_prod_fp_q1_L_4ord.dat",     branches_q1,     4},
        {"data/csv/poly_prod_q1_sv_L_4ord.dat",     branches_sv,     5} //*/ 

        // sv <= [Poly] q1-fwd <= _Poly_ <= fp 
        {"data/csv/poly_fits_fp_q1-fwd_L_4ord.dat", branches_fwd_q1, 4}, 
        {"data/csv/poly_prod_q1_sv_L_4ord.dat",     branches_sv,     5} //*/ 

        /*/ sv <= _Poly_ q1-rev <= [Poly] <= fp 
        {"data/csv/poly_prod_fp_q1_L_4ord.dat",     branches_q1,     4}, 
        {"data/csv/poly_fits_q1-rev_sv_L_4ord.dat", branches_sv,     5} //*/ 
    }; 

    const vector<NPolyArrayConstructor_t> path_and_coords_fwd{

        // sv => [Poly] => fp
        {"data/csv/poly_prod_sv_fp_L_4ord.dat", branches_fp, 5} //*/ 

        /*/ sv => _Poly_ => fp
        {"data/csv/poly_fits_sv_fp_L_4ord.dat", branches_fp, 5} //*/ 

        /*/ sv => [Poly] q1-fwd => _Poly_ => fp 
        {"data/csv/poly_prod_sv_q1_L_3ord.dat",     branches_q1, 5},
        {"data/csv/poly_fits_q1-fwd_fp_L_2ord.dat", branches_fp, 5} //*/ 

        /*/ sv <= _Poly_ q1-rev <= [Poly] <= fp 
        {"data/csv/poly_prod_fp_q1_L_4ord.dat",     branches_q1,     4}, 
        {"data/csv/poly_fits_q1-rev_sv_L_4ord.dat", branches_sv,     5} //*/ 
    };


    NPolyArrayChain chain_rev, chain_fwd; 
    
    //try to parse all polys from files
    try {
        
        for (const auto& path_and_coord : path_and_coords_rev) {

            const char* path  = path_and_coord.path.c_str(); 
            const auto  coord = path_and_coord.coords;
            const int   DoF   = path_and_coord.input_DoF;  

            chain_rev.AppendArray( ApexOptics::Parse_NPolyArray_from_file(path, coord, DoF) ); 
        }

        for (const auto& path_and_coord : path_and_coords_fwd) {

            const char* path  = path_and_coord.path.c_str(); 
            const auto  coord = path_and_coord.coords;
            const int   DoF   = path_and_coord.input_DoF;  

            chain_fwd.AppendArray( ApexOptics::Parse_NPolyArray_from_file(path, coord, DoF) ); 
        }

    } catch (const std::exception& e) {

        Error(here, "Something went wrong parsing one of the NPolyArrays.\n what(): %s", e.what()); 
        return -1; 
    }


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

    //reconstruct 'Xsv' using the 'fp=>sv' chain 
    rna.Define("Xsv_first_guess",  [&chain_rev](const Trajectory_t& Xfp)
        {
            RVec<double>&& Xfp_rvec = ApexOptics::Trajectory_t_to_RVec(Xfp); 

            RVec<double>&& Xsv_rvec = chain_rev.Eval(Xfp_rvec); 

            return ApexOptics::RVec_to_Trajectory_t(Xsv_rvec); 
        }, {"Xfp"});

    //now, use the 'sv=>fp' model to determine where the 'wiggle room' is in dp/p 
    const double d_dpp = 1e-3; 
    rna.Define("dXsv", [&chain_fwd, d_dpp](const Trajectory_t& Xsv)
        {
            //the jacobian of the fwd-model
            RMatrix&& J = chain_fwd.Jacobian( ApexOptics::Trajectory_t_to_RVec(Xsv) ); 

            RMatrix Ji(4,4, {
                J.get(0,0), J.get(0,1), J.get(0,2), J.get(0,3), 
                J.get(1,0), J.get(1,1), J.get(1,2), J.get(1,3), 
                J.get(2,0), J.get(2,1), J.get(2,2), J.get(2,3), 
                J.get(3,0), J.get(3,1), J.get(3,2), J.get(3,3)
            }); 

            Ji.Set_report_singular(false);
            
            RVec<double> J4{ 
                J.get(0,4), 
                J.get(1,4), 
                J.get(2,4), 
                J.get(3,4) 
            }; 

            RVec<double>&& dX = Ji.Solve( J4*(-1.) ); 
            
            //check for NaN / invalid result 
            if (dX.size() != 4)                return Trajectory_t{numeric_limits<double>::quiet_NaN()}; 
            for (double& x : dX) { if (x != x) return Trajectory_t{numeric_limits<double>::quiet_NaN()}; } 

            dX.push_back( 1. ); 

            return ApexOptics::RVec_to_Trajectory_t( dX ); 

        }, {"Xsv_first_guess"}); 

    rna.Define("dp", [is_RHRS, d_dpp](Trajectory_t Xsv, Trajectory_t dXsv, TVector3 vtx) 
        {
            //first, let's reconstruct the vectors in the HCS
            Trajectory_t Xhcs    = ApexOptics::SCS_to_HCS(is_RHRS, Xsv); 
            Trajectory_t Xhcs_dp = ApexOptics::SCS_to_HCS(is_RHRS, Xsv + dXsv); 

            double y0 = Xhcs.y    + (Xhcs.dydz    * vtx.z()); 
            double y1 = Xhcs_dp.y + (Xhcs_dp.dydz * vtx.z()); 
            
            return ( vtx.y() - y0 )/( y1 - y0 ); 

        }, {"Xsv_first_guess", "dXsv", "position_vtx"});
    
    rna.Define("Xsv_reco", [](Trajectory_t Xsv, Trajectory_t dXsv, double dp)
        {
            return Xsv + ( dXsv * dp ); 

        }, {"Xsv_first_guess", "dXsv", "dp"}); 

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
    
    rna = add_branch_from_Trajectory_t(rna.Get(), "Xsv_first_guess", {
        {"fg_x_sv",    &Trajectory_t::x},
        {"fg_y_sv",    &Trajectory_t::y},
        {"fg_dxdz_sv", &Trajectory_t::dxdz},
        {"fg_dydz_sv", &Trajectory_t::dydz}
    }); 

    rna = add_branch_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        {"reco_x_sv",    &Trajectory_t::x},
        {"reco_y_sv",    &Trajectory_t::y},
        {"reco_dxdz_sv", &Trajectory_t::dxdz},
        {"reco_dydz_sv", &Trajectory_t::dydz}
    }); 

    
    const int center_row = 8; 
    const int n_side_rows = 4; 

    const int center_col = 5; 
    const int n_side_cols = 4; 

    auto fit_result_ret = TestAngleReco(is_RHRS, rna.Get(), target, center_row, n_side_rows, center_col, n_side_cols); 

    AngleFitResult_t fit; 

    if (fit_result_ret) {
        fit = fit_result_ret.value(); 
    } else {
        Error(here, "Something went wrong with the fitresult."); 
        return -1; 
    }

    printf("position error: %.4e\n", fit.sigma_dydz_position );
    printf("smearing error: %.4e\n", fit.sigma_dydz_smearing );  

    printf("total error: %.4e\n", fit.sigma_dydz_overall); 

    return 0; 


    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

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
    gStyle->SetOptStat(0); 




    auto c = new TCanvas("c1", c_title, 1200, 600); 
    c->SetLeftMargin(0.12); c->SetRightMargin(0.05); 
    c->Divide(2,1, 0.01,0.01); 

    c->cd(1); hist_fg_xy->DrawCopy("col2"); 
    c->cd(2); hist_xy->DrawCopy("col2");
    
    auto c1 = new TCanvas("c2", c_title, 1200, 600); 
    c1->Divide(2,1, 0.01,0.01); 

    c1->cd(1); hist_fg_angles->DrawCopy("col2"); 
    c1->cd(2); hist_angles->DrawCopy("col2"); 

    return 0; 
}