
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

#if 0 
//A small helper class which is meant to avoid some of the awkward syntax assocaited with RDataFrame creation. 
class RDFNodeAccumulator {
private: 
    std::vector<ROOT::RDF::RNode> fNodes; 
public: 
    RDFNodeAccumulator(ROOT::RDF::RNode start) : fNodes{start} {}; 
    ~RDFNodeAccumulator() {}; 

    template<typename F> void Define(const char* new_branch, F expression, const vector<string>& inputs) {
        try {        
            fNodes.emplace_back( fNodes.back().Define(new_branch, expression, inputs) ); 
      
        } catch (const std::exception& e) {
            Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what()); 
            fNodes.clear();
            std::abort(); 
        }
    }; 

    //same as above, but this allows passing the new branch name as a std::string
    template<typename F> void Define(string new_branch, F expression, const vector<string>& inputs) {
        Define(new_branch.c_str(), expression, inputs); 
    }

    //same as above, but only define if this column is not yet defined in the dataframe
    template<typename F> void DefineIfMissing(string new_branch, F expression, const vector<string>& inputs) {

        if (!IsBranchDefined(new_branch)) Define(new_branch, expression, inputs); 
    }


    //check if a branch is defined in the dataframe as it currently exists
    bool IsBranchDefined(string branch) {
        for (const string& column : fNodes.back().GetColumnNames()) { if (branch == column) return true; }
        return false; 
    }
    //bool IsBranchDefined(const char* branch) const { string str(branch); return IsBranchDefined(str); }

    ROOT::RDF::RNode& Get() { return fNodes.back(); }

    //assignment operator 
    ROOT::RDF::RNode& operator=(ROOT::RDF::RNode node) {
        fNodes.emplace_back(node); 
        return fNodes.back(); 
    }
};
#endif 

//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int test_forward_model( const char* path_infile="data/replay/real_L_V2_sieve.root",
                        const char* path_dbfile="data/csv/poly_WireAndFoil_fp_sv_L_4ord.dat",  
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
        "dydz_sv",
        "dpp_sv"
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

    //this is a little helper class which is meant to avoid some of the awkward syntax typically associated with RDataFrames creation. 
    auto rna = RDFNodeAccumulator(df); 

    //probably not the most elegant way to do this, but here we are. 
    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{ x, y, dxdz, dydz };
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Define("Xsv_reco",  [&parr_forward](const Trajectory_t& Xfp)
        {
            return ApexOptics::RVec_to_Trajectory_t(
                parr_forward.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz})
            ); 
        }, {"Xfp"});

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
    
    rna.Overwrite("position_vtx_scs",  [&HCS_to_SCS](TVector3 vtx_hcs)
        {
            TVector3 vtx_scs = vtx_hcs; 
            HCS_to_SCS(&vtx_scs); 
            return vtx_scs; 
        }, {"position_vtx"}); 

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

    rna.Overwrite("reco_x_sv", [](double x, const TVector3& vtx){ return x + vtx.x(); }, {"reco_x_sv", "position_vtx_scs"});
    
    rna.Overwrite("reco_dxdz_sv", [](double dxdz, const TVector3& vtx){ return dxdz + (vtx.x() /( 0. - vtx.z())); }, {"reco_dxdz_sv", "position_vtx_scs"}); 

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    //create both histograms
    auto hist_xy = rna.Get()
        .Histo2D({"h_xy",     "Sieve-plane projection;x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "reco_x_sv", "reco_y_sv");
    
    auto hist_angles = rna.Get()
        .Histo2D({"h_angles", "Sieve-plane projection;dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "reco_dxdz_sv", "reco_dydz_sv"); 
    
    auto hist_z_y = rna.Get()
        .Histo2D({"h_z_y",    "z_{tg} vs y_{sv};z_{tg};y_{sv}", 200, -0.4, 0.4, 200, -0.035, 0.035}, "z_reco_vertical", "reco_y_sv");

    auto hist_z_dydz = rna.Get()
        .Histo2D({"h_z_dydz", "z_{tg} vs dy/dz_{sv};z_{tg};dy/dz_{sv}", 200, -0.4, 0.4, 200, -0.035, 0.035}, "z_reco_vertical", "reco_dydz_sv");

    auto hist_dydz = rna.Get()
        .Filter([](double dxdz){ return (0.0e-3 < dxdz) && (1.0e-3 > dxdz); }, {"reco_dxdz_sv"})
        .Histo1D<double>({"h_dydz", "dy/dz_{tg} (horizontal angle);dy/dz_{tg} (rad)", 200, -0.04, 0.03}, "reco_dydz_sv"); 


    auto hist_z_v = rna.Get()
        .Histo1D<double>({"h_z_vert", "z_{tg} vertical plane;z_{tg};", 200, -0.4, 0.4}, "z_reco_vertical");

    auto hist_z_h = rna.Get()
        .Histo1D<double>({"h_z_horiz", "z_{tg} horizontal plane;z_{tg};", 200, -0.4, 0.4}, "z_reco_horizontal");

    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 

    //PolynomialCut::InteractiveApp((TH2*)hist_z_y->Clone("hclone"), "col2", kSunset); 
    //return 0; 
    
    auto c = new TCanvas("c1", c_title, 1200, 600); 
    c->Divide(2,1, 0.01,0.01); 

    c->cd(1); hist_xy->DrawCopy("col2"); 
    c->cd(2); hist_angles->DrawCopy("col2");
    
    auto c1 = new TCanvas("c2", c_title, 1200, 600); 
    c1->Divide(2,1, 0.01,0.01); 

    c1->cd(1); hist_z_y->DrawCopy("col2"); 
    c1->cd(2); hist_z_dydz->DrawCopy("col2"); 

    auto c2 = new TCanvas("c3", c_title); 
    hist_dydz->DrawCopy(); 

    return 0; 
}