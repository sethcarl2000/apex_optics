
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include <vector>
#include <array> 
#include <chrono> 
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

using namespace std; 
using namespace ROOT::VecOps; 


struct Track_t { 
    double x,y,dxdz,dydz,dpp; 
}; 

std::vector<double> *pt_x, *pt_y; 
bool *clicked; 

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
int fwdmodel_pick_holes( const char* path_infile="data/replay/real_L_V2_sieve.root",
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
    
    array<double,2> hole_center{+2.22055e-02, -3.17677e-02}; 
    double radius = 3e-3; 

    
    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 


    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //probably not the most elegant way to do this, but here we are. 
    
    auto df_reco = df 

        .Define("position_vtx_scs", [&HCS_to_SCS](TVector3 vtx_hcs){
            TVector3 vtx_scs = vtx_hcs; 
            HCS_to_SCS(&vtx_scs); 
            return vtx_scs; 
        }, {"position_vtx"})

        .Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Track_t{ .x=x, .y=y, .dxdz=dxdz, .dydz=dydz }; 
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"})

        .Define("Xsv",  [&parr_forward](const Track_t& Xfp)
        {
            auto Xsv = parr_forward.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz});
            return Track_t{ .x=Xsv[0], .y=Xsv[1], .dxdz=Xsv[2], .dydz=Xsv[3] }; 
        }, {"Xfp"})

        .Define("Xsv_rvtx", [](Track_t Xsv, TVector3 vtx_scs)
        {
            Xsv.x += vtx_scs.x(); 
            Xsv.y += vtx_scs.y(); 
            return Xsv; 
        }, {"Xsv", "position_vtx_scs"})

        .Filter([hole_center, radius](const Track_t& Xsv)
        {
            return pow(Xsv.x - hole_center[0], 2) + pow(Xsv.y - hole_center[1], 2) < radius*radius; 
        }, {"Xsv_rvtx"}); 

        

    auto df_out = add_branch_from_Track_t(df_reco, "Xsv", {
        {"x_sv",    &Track_t::x},
        {"y_sv",    &Track_t::y},
        {"dxdz_sv", &Track_t::dxdz},
        {"dydz_sv", &Track_t::dydz}
    }); 

    //create both histograms
    auto hist_xy = df_out
        
        //correct for react-vertex position
        .Define("x_react_vtx_fix", [](double x, TVector3 r){return x + r.x();}, {"x_sv", "position_vtx_scs"})
        .Define("y_react_vtx_fix", [](double y, TVector3 r){return y + r.y();}, {"y_sv", "position_vtx_scs"})

        .Histo2D({"h_xy", "Sieve-plane projection;x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, 
                "x_react_vtx_fix", "y_react_vtx_fix");
    
    auto hist_angles    
        = df_out.Histo2D({"h_angles", "Sieve-plane projection;dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "dxdz_sv", "dydz_sv"); 

    
    vector<string> branches_input{
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    vector<string> branches_output{
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv"
    }; 
    
    const size_t n_input = branches_input.size(); 
    
    
    
    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 
    
    vector<TCanvas*> v_canv; 
    
    int i_page=0; 
    for (const string& br_input : branches_input) {
        
        //NOTE: assumes that the canvas has 4 pannels!!!!!
        int i_canv=1;

        char c_name[50]; sprintf(c_name, "c_%i", i_page); 

        v_canv.push_back(new TCanvas(c_name, c_title, 1200, 800)); 
        auto &canv = v_canv.back(); 
        canv->Divide(2,2, 0.01,0.01); 

        for (const string& br_output : branches_output) {

            char h_name[200];   sprintf(h_name, "h_%s_%s", br_output.c_str(), br_input.c_str()); 
            char h_title[200];  sprintf(h_title, ";%s;%s", br_input.c_str(), br_output.c_str()); 
            
            auto x_min = df_out.Min(br_input); 
            auto x_max = df_out.Max(br_input); 
            
            auto y_min = df_out.Min(br_output); 
            auto y_max = df_out.Max(br_output); 
             
            auto hist = df_out.Histo2D<double>({h_name, h_title, 200, *x_min, *x_max, 200, *y_min, *y_max}, br_input, br_output); 
            
            canv->cd(i_canv++); 
            hist->DrawCopy("col2"); 

            canv->Modified(); 
            canv->Update(); 
        }
        
        printf("done with page %i \n.", i_page++); 

        /*canv->cd(0); 
        switch (i_page) {
            case 1:       canv->SaveAs("correlations.pdf("); break;  //open the pdf
            case 4:       canv->SaveAs("correlations.pdf)"); break;  //close the pdf 
            default:      canv->SaveAs("correlations.pdf");          //intermediate page
        } */ 
        
    }


    
    //now, start the event loop; 
    pt_x = new vector<double>{};
    pt_y = new vector<double>{}; 
    
    //clicked = new bool(false); 
    //gPad->AddExec("x1", ".x scripts/select_hole.C(\"test\")"); 

    return 0; 
}