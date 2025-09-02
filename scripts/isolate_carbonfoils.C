
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TVector3.h> 
#include <ApexOptics.h> 
#include <TSystem.h>
#include <TParameter.h> 
#include <PolynomialCut.h> 
#include <TStyle.h>
#include <TCanvas.h> 
#include "include/RDFNodeAccumulator.h"

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
int isolate_carbonfoils(    const char* path_infile ="data/replay/real_L_Opt1-4773.root",
                            const char* target_name="none",       
                            const char* path_cutfile="data/csv/polycut_L_O1-new.dat",
                            const char* path_outfile="data/replay/real_L_O1.root",
                            const char* path_dbfile ="data/csv/poly_WireAndFoil_fp_sv_L_4ord.dat",  
                            const char* tree_name   ="tracks_fp" ) 
{
    const char* const here = "test_forward_model"; 

    auto infile = new TFile(path_infile, "READ");
    
    OpticsTarget_t target; 
    try { target = ApexOptics::GetTarget(string(target_name)); } 

    catch (const std::exception& e) {

        Error(here, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
        return -1; 
    }
    
    printf("Optics target chosen is '%s'.\n", target.name.c_str()); 

    
    const bool create_new_cut = (target.name == "none"); 

    if (create_new_cut) { cout << "Creating a new polycut." << endl; }


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

    //try to create the polynomial cut
    auto polycut = new PolynomialCut; 

    if (!create_new_cut) {
        try { 
            polycut->Parse_dbfile(path_cutfile); 
        } 
        catch (const PolynomialCut::DBFileException& e) {
            Error(here, "Problem parsing cutfile. Exception:\n %s", e.what());
            return -1;  
        }
    }
    const char* cutbranch_x = "z_reco_vertical"; 
    const char* cutbranch_y = "y_sv"; 
    
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

    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 


    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 
    
    //create the 'RDFNodeAccumulator' object. 
    RDFNodeAccumulator rna(df); 

    //probably not the most elegant way to do this, but here we are. 
    
    rna.Define("Xfp",       [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{ .x=x, .y=y, .dxdz=dxdz, .dydz=dydz }; 
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"}); 

    rna.Define("Xsv_reco",  [&parr_forward](const Trajectory_t& Xfp)
        {
            auto Xsv = parr_forward.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz});
            return Trajectory_t{ .x=Xsv[0], .y=Xsv[1], .dxdz=Xsv[2], .dydz=Xsv[3] }; 
        }, {"Xfp"}); 

    rna.Define("Xhcs_reco", [is_RHRS](const Trajectory_t& Xsv)
        {
            return ApexOptics::SCS_to_HCS(is_RHRS, Xsv); 
        }, {"Xsv_reco"});

    rna.Define("z_reco_horizontal", [](const Trajectory_t& Xhcs, const TVector3& vtx)
        {   
            return - ( Xhcs.y - vtx.y() ) / Xhcs.dydz; 
        }, {"Xhcs_reco", "position_vtx"});

    rna.Define("z_reco_vertical",   [](const Trajectory_t& Xhcs, const TVector3& vtx)
        {
            return - ( Xhcs.x - vtx.x() ) / Xhcs.dxdz; 
        }, {"Xhcs_reco", "position_vtx"});

    //this is just to keep the 'old' value, to see how good our cut is. 
    rna.Define("y_vtx", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}); 

    rna.Overwrite("position_vtx", [&target](TVector3 vtx)
        {
            //overwrite whichever coordinates are precisely given by this target's geometry. 
            //if any of the target's coordinates are NAN, then they will NOT overwrite that coordinate of the react-vetex. 
            if (target.x_hcs == target.x_hcs) vtx[0] = target.x_hcs;
            if (target.y_hcs == target.y_hcs) vtx[1] = target.y_hcs;
            if (target.z_hcs == target.z_hcs) vtx[2] = target.z_hcs;

            return vtx; 
        }, {"position_vtx"}); 
        
    rna.Define("position_vtx_scs", [is_RHRS](TVector3 vtx)
        {
            return ApexOptics::HCS_to_SCS(is_RHRS, vtx);
        }, {"position_vtx"});


    rna = add_branch_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        {"x_sv",    &Trajectory_t::x},
        {"y_sv",    &Trajectory_t::y},
        {"dxdz_sv", &Trajectory_t::dxdz},
        {"dydz_sv", &Trajectory_t::dydz},
        {"dpp_sv",  &Trajectory_t::dpp}
    });  

    auto node_precut = rna.Get(); 

    if (!create_new_cut) { 
        rna = rna.Get()
            .Filter([polycut](double cut_x, double cut_y)
            { 
                return polycut->IsInside(cut_x, cut_y); 
            }, {cutbranch_x, cutbranch_y}); 
    }

    //if we aren't creating a new cut, then let's make a snapshot (create a new ROOT file)
    if (!create_new_cut) {
        
        rna.Get()
            .Filter([polycut](double cut_x, double cut_y)
            { 
                return polycut->IsInside(cut_x, cut_y); 
            }, {cutbranch_x, cutbranch_y})

            .Snapshot("tracks_fp", path_outfile, {
                "x_sv",
                "y_sv",
                "dxdz_sv",
                "dydz_sv",
                "dpp_sv",

                "x_fp",
                "y_fp",
                "dxdz_fp",
                "dydz_fp",

                "position_vtx", 
                "position_vtx_scs"
            });              
    }


    //create both histograms
    auto hist_xy = rna.Get()
        //correct for react-vertex position
        /*.Define("x_react_vtx_fix", [](double x, TVector3 r){return x + r.x();}, {"x_sv", "position_vtx_scs"})
        .Define("y_react_vtx_fix", [](double y, TVector3 r){return y + r.y();}, {"y_sv", "position_vtx_scs"})*/ 

        .Histo2D({"h_xy", "Sieve-plane projection;x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "x_sv", "y_sv");
    
    auto hist_angles = rna.Get()
        .Histo2D({"h_angles", "Sieve-plane projection;dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "dxdz_sv", "dydz_sv"); 
    
    
    auto hist_z_y = rna.Get()
        .Histo2D({"h_z_y", "z_{tg} vs y_{sv};z_{tg};y_{sv}", 200, -0.45, 0.45, 200, -0.045, 0.045}, "z_reco_vertical", "y_sv");

    auto hist_z_dydz = rna.Get()
        .Histo2D({"h_z_dydz", "z_{tg} vs dy/dz_{sv};z_{tg};dy/dz_{sv}", 200, -0.4, 0.4, 200, -0.04, 0.04}, "z_reco_vertical", "dydz_sv");


    auto hist_y_vtx_precut = node_precut 
        .Histo1D({"h_all", "All events;y_{hcs} (m);", 200, -3e-3, 8e-3}, "y_vtx"); 

    auto hist_y_vtx_postcut = rna.Get()
        .Histo1D({"h_cut", "in polynomial cut;y_{hcs} (m);", 200, -3e-3, 8e-3}, "y_vtx"); 

    char c_title[256]; 
    sprintf(c_title, "data:'%s', cutfile:'%s'", path_infile, path_cutfile); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 

    if (create_new_cut) {
        
        PolynomialCut::InteractiveApp((TH2*)hist_z_y->Clone("hclone"), "col2", kSunset); 
    
    } else {



        auto c = new TCanvas("c1", c_title, 1200, 600); 
        c->Divide(2,1, 0.01,0.01); 

        c->cd(1); hist_xy->DrawCopy("col2"); 
        c->cd(2); hist_angles->DrawCopy("col2");
        
        auto c1 = new TCanvas("c2", c_title, 1200, 600); 
        c1->Divide(2,1, 0.01,0.01); 

        c1->cd(1); hist_z_y->DrawCopy("col2"); 
        c1->cd(2); hist_z_dydz->DrawCopy("col2"); 

        
        auto c2 = new TCanvas("c3", c_title); 
        hist_y_vtx_precut->DrawCopy(); 
        
        hist_y_vtx_postcut->SetFillStyle(3004);
        hist_y_vtx_postcut->SetFillColor(kRed); 
        hist_y_vtx_postcut->SetLineColor(kRed); 

        hist_y_vtx_postcut->DrawCopy("SAME"); 

        gPad->BuildLegend(); 




            
        //write 'is_RHRS' parameter 
        auto file = new TFile(path_outfile, "UPDATE"); 
        param_is_RHRS = new TParameter<bool>("is_RHRS", is_RHRS); 
        
        param_is_RHRS->Write(); 
        file->Close();
    }
    return 0; 
}