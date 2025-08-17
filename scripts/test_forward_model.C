
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

using namespace std; 
using namespace ROOT::VecOps; 


struct Track_t { 
    double x,y,dxdz,dydz,dpp; 
}; 

NPoly dummy; 

//transform coordinates from Sieve Coordinate System (SCS) to Hall coordinate system (HCS)
void SCS_to_HCS(Track_t& track, const bool is_RHRS) 
{
    //direction (SCS)
    auto dir = TVector3( track.dxdz, track.dydz, 1. );

    auto pos = TVector3( track.x, track.y, 0. ) + ApexOptics::Get_sieve_pos(is_RHRS); 

    //rotate both the position and the direction
    dir.RotateZ( -TMath::Pi()/2. ); 
    dir.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 

    pos.RotateZ( -TMath::Pi()/2. ); 
    pos.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 

    //compute the new slopes
    track.dxdz = dir.x() / dir.z(); 
    track.dydz = dir.y() / dir.z(); 

    //use these new slopes to project the track onto the z=0 plane in HCS 
    track.x = pos.x() - track.dxdz * pos.z(); 
    track.y = pos.y() - track.dydz * pos.z(); 
}

//_______________________________________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Track_t struct. 
ROOT::RDF::RNode add_branch_from_Track_t(ROOT::RDF::RNode df, const char* branch_in, map<string, double Track_t::*> branches)
{
    const int n_nodes = branches.size() + 1; 
    RVec<ROOT::RDF::RNode> nodes{ df }; 

    int i_branch=0; 
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double Track_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Track_t::coord' 
        nodes.push_back( nodes.back()
            
            .Define(branch_name, [coord](const Track_t& track) { return track.*coord; }, {branch_in})
        ); 
    }
    
    return nodes.back(); 
}



//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int test_forward_model( const char* path_infile="data/replay/real_L_V2_sieve.root",
                        const char* path_dbfile="data/csv/db_real_fp_sv_v2_L_2ord.dat",  
                        const char* tree_name="tracks_fp" ) 
{
    const char* const here = "test_forward_model"; 

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

    
    vector<string> branches_sv{
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv"
    }; 

    NPolyArray parr_forward = ApexOptics::Parse_NPolyArray_from_file(path_dbfile, branches_sv, 4); 
    
    //check to see if the poly was parsed successfully
    if (parr_forward.Get_status() != NPolyArray::kGood) {
        Error(here, "Problem parsing NPolyArray from file '%s'", path_dbfile); 
        return -1; 
    }

    //transform coordinates from Sieve Coordinate System (SCS) to Hall coordinate system (HCS)
    auto HCS_to_SCS = [is_RHRS](TVector3 *pos, TVector3 *dir=nullptr) 
    {
        //rotate both the position and the direction
        /*
        dir.RotateZ( -TMath::Pi()/2. ); 
        dir.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 

        pos.RotateZ( -TMath::Pi()/2. ); 
        pos.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 
        */ 

        pos->RotateY( -ApexOptics::Get_sieve_angle(is_RHRS) ); 
        pos->RotateZ( TMath::Pi()/2. ); 

        *pos += - ApexOptics::Get_sieve_pos(is_RHRS);
        
        if (dir) {
            dir->RotateY( -ApexOptics::Get_sieve_angle(is_RHRS) ); 
            dir->RotateZ( TMath::Pi()/2. ); 
        }
        return; 
    }; 

    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 


    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //probably not the most elegant way to do this, but here we are. 
    
    auto df_reco = df 

        .Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Track_t{ .x=x, .y=y, .dxdz=dxdz, .dydz=dydz }; 
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"})

        .Define("Xsv_reco",  [&parr_forward](const Track_t& Xfp)
        {
            auto Xsv = parr_forward.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz});
            return Track_t{ .x=Xsv[0], .y=Xsv[1], .dxdz=Xsv[2], .dydz=Xsv[3] }; 
        }, {"Xfp"})

        .Define("Xhcs_reco", [is_RHRS](const Track_t& Xsv)
        {
            Track_t Xhcs{Xsv}; 
            SCS_to_HCS(Xhcs, is_RHRS); 
            return Xhcs; 
        }, {"Xsv_reco"})

        .Define("z_reco_horizontal", [](const Track_t& Xhcs, const TVector3& vtx)
        {   
            return - ( Xhcs.y - vtx.y() ) / Xhcs.dydz; 
        }, {"Xhcs_reco", "position_vtx"})

        .Define("z_reco_vertical",   [](const Track_t& Xhcs, const TVector3& vtx)
        {
            return - ( Xhcs.x - vtx.x() ) / Xhcs.dxdz; 
        }, {"Xhcs_reco", "position_vtx"})

        .Define("position_vtx_scs", [&HCS_to_SCS](TVector3 vtx_hcs){
            TVector3 vtx_scs = vtx_hcs; 
            HCS_to_SCS(&vtx_scs); 
            return vtx_scs; 
        }, {"position_vtx"}); 

    auto df_out = add_branch_from_Track_t(df_reco, "Xsv_reco", {
        {"reco_x_sv",    &Track_t::x},
        {"reco_y_sv",    &Track_t::y},
        {"reco_dxdz_sv", &Track_t::dxdz},
        {"reco_dydz_sv", &Track_t::dydz}
    }); 

    //create both histograms
    auto hist_xy = df_out
        
        //correct for react-vertex position
        .Define("x_react_vtx_fix", [](double x, TVector3 r){return x + r.x();}, {"reco_x_sv", "position_vtx_scs"})
        .Define("y_react_vtx_fix", [](double y, TVector3 r){return y + r.y();}, {"reco_y_sv", "position_vtx_scs"})

        .Histo2D({"h_xy", "Sieve-plane projection;x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, 
                "x_react_vtx_fix", "y_react_vtx_fix");
    
    auto hist_angles    
        = df_out.Histo2D({"h_angles", "Sieve-plane projection;dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "reco_dxdz_sv", "reco_dydz_sv"); 
    
    
    auto hist_y_dydz 
        = df_out.Histo2D({"h_y_dydz", "y_{sv} vs dy/dz_{sv};y_{sv};dy/dz_{sv}", 200, -0.07, 0.07, 200, -0.035, 0.035}, "reco_y_sv", "reco_dydz_sv");

    auto hist_z_dydz 
        = df_out.Histo2D({"h_x_dxdz", "z_{tg} vs y_{sv};z_{tg};y_{sv}", 200, -0.4, 0.4, 200, -0.035, 0.035}, "z_reco_vertical", "reco_y_sv");


    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 
    
    auto c = new TCanvas("c1", c_title, 1200, 600); 
    c->Divide(2,1, 0.01,0.01); 

    c->cd(1); hist_xy->DrawCopy("col2"); 
    c->cd(2); hist_angles->DrawCopy("col2");
    
    auto c1 = new TCanvas("c2", c_title, 1200, 600); 
    c1->Divide(2,1, 0.01,0.01); 

    c1->cd(1); hist_y_dydz->DrawCopy("col2"); 
    c1->cd(2); hist_z_dydz->DrawCopy("col2"); 

    return 0; 
}