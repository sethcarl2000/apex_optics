
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
//this last one should eventually be moved to the ApexOptics namespace, as it has general use outside of the 'isolate_sieveholes' app. 
#include <isolate_sieveholes/SieveHoleData.h>
#include <TLine.h> 

using namespace std; 
using namespace ROOT::VecOps; 

using ApexOptics::Trajectory_t; 
using ApexOptics::OpticsTarget_t; 

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


//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int test_angular_accuracy(  const char* path_infile ="data/replay/real_L_V2.root",
                            std::string target_name ="V2",
                            const char* path_dbfile ="data/csv/poly_WireAndFoil_fp_sv_L_4ord.dat",  
                            const char* path_dXsv   ="data/csv/poly_dXsv_L_2ord.dat",
                            const char* tree_name   ="tracks_fp" ) 
{
    const char* const here = "test_angular_accuracy"; 

    //try to get the target we need. 
    OpticsTarget_t target; 
    try { 
        target = ApexOptics::GetTarget(string(target_name)); 
    
    } catch (const std::exception& e) {

        Error(here, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
        return -1; 
    }

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

    //are we testing the vertical or horizontal angle? 
    // if x_hcs is NOT nan, that means that this is a vertical wire run. 
    // otherwise, its a horizontal-wire run. 
    const bool test_horizontal_angle = (target.x_hcs == target.x_hcs); 


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

    //choose which row/column to start at, and how many rows/columns to do in either direction. 

    const int center_row = 8; 
    const int n_side_rows = 3; 

    const int center_col = 7; 
    const int n_side_cols = 3; 



    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 

    

    const double y_rast_rms = 0.113e-3; 

    const double d_dpp = 1e-3; 

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //get the average react-vertex. 
    TVector3 react_vtx_hcs(
        0.,
        0.,
        target.z_hcs
    ); 

    //figure out which coordinates are well-defined for this target (they are 'nan' if not). 
    if (target.x_hcs == target.x_hcs) {
        react_vtx_hcs[0] = target.x_hcs; 
    } else {
        react_vtx_hcs[0] = *df.Define("x", [](TVector3 vtx){ return vtx.x(); }, {"position_vtx"}).Mean("x"); 
    }

    if (target.y_hcs == target.y_hcs) {
        react_vtx_hcs[1] = target.y_hcs; 
    } else {
        react_vtx_hcs[1] = *df.Define("y", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}).Mean("y"); 
    }
    
    //get the react vertex in Sieve Coordinates
    const TVector3 react_vtx_scs = ApexOptics::HCS_to_SCS(is_RHRS, react_vtx_hcs); 

    



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
            return ApexOptics::RVec_to_Trajectory_t( parr_fp_sv.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz}) ); 

        }, {"Xfp"});

    //define a y-react pos. if it's monte-carlo data, then we want to add a realistic RMS. 
    rna.Define("y_vtx", [y_rast_rms, &rand, is_MonteCarlo](TVector3 vtx)
        {
            return vtx.y() + ( is_MonteCarlo ? rand.Gaus() * y_rast_rms : 0. ); 

        }, {"position_vtx"}); 

    rna.Define("dXsv",     [&parr_dXsv](Trajectory_t Xsv)
        {
            return ApexOptics::RVec_to_Trajectory_t( parr_dXsv.Eval({Xsv.x, Xsv.y, Xsv.dxdz, Xsv.dydz, Xsv.dpp}) ); 
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
            return Xsv + ( dXsv * dp ); 

        }, {"Xsv_first_guess", "dXsv", "dp"});

    rna.Define("Xsv_reco_upfoil", [](Trajectory_t Xsv, Trajectory_t dXsv, double dp)
        {
            return Xsv + ( dXsv * dp ); 

        }, {"Xsv_first_guess", "dXsv", "dp_upfoil"});

    rna.Define("Xsv_reco_downfoil", [](Trajectory_t Xsv, Trajectory_t dXsv, double dp)
        {
            return Xsv + ( dXsv * dp ); 

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


    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 


    //now, we can compute where each sieve-hole SHOULD be
    vector<SieveHole> sieve_holes = SieveHole::ConstructSieveHoles(is_RHRS); 

    auto FindSieveHole = [&sieve_holes](int row, int col) 
    {
        auto it = std::find_if( 
            sieve_holes.begin(),
            sieve_holes.end(), 
            [row,col](const SieveHole& elem){ return row==elem.row && col==elem.col; }
        ); 
        if (it == sieve_holes.end()) {
            ostringstream oss; 
            oss << "in <FindSieveHole>: invalid sieve hole requested. row/col: " << row << "/" << col; 
            throw invalid_argument(oss.str()); 
            return SieveHole{}; 
        } 
        return SieveHole{*it}; 
    }; 

    auto c = new TCanvas("c", c_title, 1600, 800); 
    c->Divide( 2,1, 0.01,0.01 ); 

    const int canvId_x_y = 1;
    const int canvId_dx_dy = 2; 

    gStyle->SetPalette(kBird);

    if (test_horizontal_angle) {
        //here, we're testing vertical wires. 
        
        //c->Divide( 1 + n_side_rows*2, 1, 0.01,0.); 
        int i_canv=1; 
        

        auto h_x_y = rna.Get() 
            .Histo2D<double>({"h", "x_{sv} vs y_{sv} (m);", 200, -0.04, 0.04, 200, -0.04, 0.04}, "reco_x_sv", "reco_y_sv"); 

        c->cd(canvId_x_y); 
        h_x_y->DrawCopy("col2"); 


        auto h_dx_dy = rna.Get()
            .Histo2D<double>({"h", "dx/dz_{sv} vs dy/dz_{sv} (rad)", 200, -0.04, 0.04, 200, -0.04, 0.04}, "reco_dxdz_sv", "reco_dydz_sv"); 

        c->cd(canvId_dx_dy); 
        h_dx_dy->DrawCopy("col2"); 


        //loop thru eaach column 
        for (int row = center_row-n_side_rows; row < center_row+n_side_rows; row++) {

            vector<SieveHole> row_holes; 
            for (int col=center_col-n_side_cols; col<center_col+n_side_cols; col++) {
                try {
                    row_holes.push_back( FindSieveHole(row, col) ); 
                } catch (const std::exception& e) {
                    Error(here, "Error trying to access sieve hole.\n what(): %s", e.what()); 
                    return -1; 
                }
            }

            vector<Trajectory_t> hole_angles;
            for (auto& hole : row_holes) { 
                hole_angles.emplace_back(Trajectory_t{
                    hole.x, 
                    hole.y,
                    ( hole.x - react_vtx_scs.x() )/( 0. - react_vtx_scs.z() ),
                    ( hole.y - react_vtx_scs.y() )/( 0. - react_vtx_scs.z() ) 
                }); 
            }
            
            //..
            const double row_x    = hole_angles[0].x; 
            
            TLine *line; 

            c->cd(canvId_x_y); 
            line = new TLine( row_x , -0.04, row_x, 0.04 ); 
            line->Draw("SAME");
            
            //..
            const double row_dxdz = hole_angles[0].dxdz; 
            
            c->cd(canvId_dx_dy); 
            line = new TLine( row_dxdz , -0.04, row_dxdz, 0.04 ); 
            line->Draw("SAME");
        }
    }


    return 0; 
}